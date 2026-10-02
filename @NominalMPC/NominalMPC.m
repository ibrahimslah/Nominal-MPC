classdef NominalMPC < handle
    %NOMINALMPC Stateful path-coordinate MPC controller for an AUV.
    %   Construct this controller after generating corridors for the same
    %   path used by pathModel. The caller owns plant stepping and logging.
    %   The corridor manager must support getPathConstraintsAtS(s, frameSegment).
    %   info.FrameSegmentIndex identifies the frame of PathState and all
    %   PredictedState samples; info.SegmentIndex identifies the corridor.

    properties (SetAccess = private)
        PathModel
        CorridorManager
        Nx
        Nu
        N
        Ts
        Q
        R
        ReferenceSpeed
        Umin
        Umax
        GoalTolerance
        SolverOptions
        PreviousInput
        PreviousSolution
    end

    methods
        function obj = NominalMPC(pathModel, corridorManager, config)
            %NOMINALMPC Set up the controller and validate its dependencies.
            if nargin ~= 3 || ~isa(pathModel, 'AUVPathModel')
                error('NominalMPC:InvalidPathModel', ...
                    'Provide an AUVPathModel and a generated corridor manager.');
            end
            if ~isobject(corridorManager) || ...
                    ~isprop(corridorManager, 'Path') || ...
                    ~isprop(corridorManager, 'Corridors') || ...
                    ~ismethod(corridorManager, 'getPathConstraintsAtS')
                error('NominalMPC:InvalidCorridorManager', ...
                    ['Corridor manager must expose Path, Corridors, and ' ...
                    'getPathConstraintsAtS(s, frameSegment).']);
            end
            if ~isequal(corridorManager.Path, pathModel.Path)
                error('NominalMPC:PathMismatch', ...
                    'Corridor manager and path model must use the same path.');
            end
            if numel(corridorManager.Corridors) ~= size(pathModel.Path, 1) - 1
                error('NominalMPC:CorridorsNotGenerated', ...
                    'Generate one corridor per path segment before constructing NominalMPC.');
            end
            if ~isstruct(config) || ~isscalar(config)
                error('NominalMPC:InvalidConfig', 'config must be a scalar struct.');
            end

            required = {'N', 'Q', 'R', 'ReferenceSpeed', 'Umin', 'Umax'};
            for i = 1:numel(required)
                if ~isfield(config, required{i})
                    error('NominalMPC:InvalidConfig', ...
                        'Missing config field %s.', required{i});
                end
            end

            obj.PathModel = pathModel;
            obj.CorridorManager = corridorManager;
            obj.Nx = pathModel.nx;
            obj.Nu = pathModel.nu;
            obj.Ts = pathModel.SampleTime;
            validateattributes(obj.Ts, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'positive'});

            validateattributes(config.N, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', 'positive'});
            obj.N = double(config.N);

            validateattributes(config.Q, {'numeric'}, ...
                {'size', [obj.Nx obj.Nx], 'real', 'finite'});
            validateattributes(config.R, {'numeric'}, ...
                {'size', [obj.Nu obj.Nu], 'real', 'finite'});
            obj.Q = double(config.Q);
            obj.R = double(config.R);

            qTolerance = 1e-10 * max(1, norm(obj.Q, 'fro'));
            rTolerance = 1e-10 * max(1, norm(obj.R, 'fro'));
            if norm(obj.Q - obj.Q.', 'fro') > qTolerance || ...
                    min(eig((obj.Q + obj.Q.') / 2)) < -qTolerance
                error('NominalMPC:InvalidConfig', ...
                    'Q must be symmetric positive semidefinite.');
            end
            [~, rCholFlag] = chol((obj.R + obj.R.') / 2);
            if norm(obj.R - obj.R.', 'fro') > rTolerance || ...
                    rCholFlag ~= 0
                error('NominalMPC:InvalidConfig', ...
                    'R must be symmetric positive definite.');
            end

            validateattributes(config.ReferenceSpeed, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'positive'});
            obj.ReferenceSpeed = double(config.ReferenceSpeed);

            validateattributes(config.Umin, {'numeric'}, ...
                {'vector', 'numel', obj.Nu, 'real', 'finite'});
            validateattributes(config.Umax, {'numeric'}, ...
                {'vector', 'numel', obj.Nu, 'real', 'finite'});
            obj.Umin = double(config.Umin(:));
            obj.Umax = double(config.Umax(:));
            if any(obj.Umin >= obj.Umax) || ...
                    any(obj.Umin > 0) || any(obj.Umax < 0)
                error('NominalMPC:InvalidConfig', ...
                    'Input bounds must be ordered and contain zero.');
            end

            obj.GoalTolerance = 0.03;
            if isfield(config, 'GoalTolerance')
                validateattributes(config.GoalTolerance, {'numeric'}, ...
                    {'scalar', 'real', 'finite', 'nonnegative'});
                obj.GoalTolerance = double(config.GoalTolerance);
            end

            if isfield(config, 'SolverOptions') && ~isempty(config.SolverOptions)
                obj.SolverOptions = config.SolverOptions;
                if ~isa(obj.SolverOptions, 'optim.options.Quadprog') || ...
                        ~strcmpi(obj.SolverOptions.Algorithm, 'active-set')
                    error('NominalMPC:InvalidConfig', ...
                        'SolverOptions must be quadprog options with Algorithm active-set.');
                end
            else
                obj.SolverOptions = optimoptions('quadprog', ...
                    'Algorithm', 'active-set', 'Display', 'off');
            end

            obj.reset();
        end

        function [uk, info] = computeControl(obj, xCartesian)
            %COMPUTECONTROL Return the first MPC input for a Cartesian state.
            validateattributes(xCartesian, {'numeric'}, ...
                {'vector', 'numel', obj.Nx, 'real', 'finite'});
            [xk, frameSegment] = obj.PathModel.updateFromCartesian(xCartesian(:));
            info = obj.emptyInfo(xk);
            % Every state in this solve uses this fixed segment frame.
            info.FrameSegmentIndex = frameSegment;

            if xk(1) >= obj.PathModel.TotalLength - obj.GoalTolerance
                uk = zeros(obj.Nu, 1);
                obj.reset();
                info.GoalReached = true;
                info.SegmentIndex = obj.PathModel.ActiveSegment;
                return
            end

            [Ad, Bd, rd] = obj.PathModel.linearizedDiscrete(xk, obj.PreviousInput);
            [Alift, Blift, rlift] = obj.liftedAffine(Ad, Bd, rd, obj.N);
            Xref = obj.buildReferenceTrajectory(xk);

            [Apath, bpath, segmentIndex] = ...
                obj.CorridorManager.getPathConstraintsAtS(xk(1), frameSegment);
            if ~isnumeric(Apath) || size(Apath, 2) ~= obj.Nx || ...
                    ~isnumeric(bpath) || ~iscolumn(bpath) || ...
                    numel(bpath) ~= size(Apath, 1) || ...
                    ~isreal(Apath) || ~isreal(bpath) || ...
                    any(~isfinite(Apath(:))) || any(~isfinite(bpath(:)))
                error('NominalMPC:InvalidCorridorConstraints', ...
                    'Corridor constraints must have dimensions m-by-Nx and m-by-1.');
            end

            [H, f] = obj.buildNominalMPCCost( ...
                xk, Alift, Blift, rlift, Xref, obj.Q, obj.R, obj.N);
            [Aineq, bineq] = obj.buildPathMPCConstraints( ...
                xk, Alift, Blift, rlift, Apath, bpath, obj.N);

            lb = repmat(obj.Umin, obj.N, 1);
            ub = repmat(obj.Umax, obj.N, 1);
            x0 = obj.initialGuess();
            [Uopt, cost, exitflag, output] = quadprog( ...
                H, f, Aineq, bineq, [], [], lb, ub, x0, obj.SolverOptions);

            info.SegmentIndex = segmentIndex;
            info.InsideCorridor = all(Apath*xk <= bpath + 1e-8);
            info.Reference = Xref;
            info.InitialGuess = x0;
            info.ExitFlag = exitflag;
            info.SolverOutput = output;
            if ~isempty(cost)
                info.Cost = cost;
            end

            if exitflag > 0
                if numel(Uopt) ~= obj.N*obj.Nu || any(~isfinite(Uopt))
                    error('NominalMPC:InvalidSolverOutput', ...
                        'quadprog reported success without a valid input sequence.');
                end
                uk = Uopt(1:obj.Nu);
                obj.PreviousSolution = obj.shiftSequence(Uopt);
                info.OptimalSequence = Uopt;
                info.PredictedState = Alift*xk + Blift*Uopt + rlift;
            else
                info.UsedFallback = true;
                if isempty(obj.PreviousSolution)
                    uk = zeros(obj.Nu, 1);
                else
                    uk = obj.PreviousSolution(1:obj.Nu);
                    obj.PreviousSolution = obj.shiftSequence(obj.PreviousSolution);
                end
            end

            obj.PreviousInput = uk;
        end

        function reset(obj)
            %RESET Clear control memory, leaving model and corridors intact.
            obj.PreviousInput = zeros(obj.Nu, 1);
            obj.PreviousSolution = [];
        end
    end

    methods (Access = private)
        [H, F] = buildNominalMPCCost(obj, xk, Alift, Blift, rlift, ...
            Xref, Q, R, N)
        [Aineq, bineq] = buildPathMPCConstraints(obj, xk, Alift, ...
            Blift, rlift, A_path, b_path, N)
        [Alift, Blift, rlift] = liftedAffine(obj, A, B, r, N)
    end

    methods (Access = private)
        function Xref = buildReferenceTrajectory(obj, xk)
            Xref = zeros((obj.N+1)*obj.Nx, 1);
            for j = 0:obj.N
                sReference = min( ...
                    xk(1) + j*obj.Ts*obj.ReferenceSpeed, ...
                    obj.PathModel.TotalLength);
                rows = j*obj.Nx + (1:obj.Nx);
                Xref(rows) = obj.PathModel.referenceState( ...
                    sReference, obj.ReferenceSpeed);
            end
        end

        function x0 = initialGuess(obj)
            if isempty(obj.PreviousSolution)
                x0 = repmat(obj.PreviousInput, obj.N, 1);
            else
                x0 = obj.PreviousSolution;
            end
        end

        function sequence = shiftSequence(obj, sequence)
            controls = reshape(sequence, obj.Nu, obj.N);
            sequence = [controls(:, 2:end), controls(:, end)];
            sequence = sequence(:);
        end

        function info = emptyInfo(~, xk)
            info = struct( ...
                'PathState', xk, ...
                'FrameSegmentIndex', [], ...
                'SegmentIndex', [], ...
                'InsideCorridor', [], ...
                'Reference', [], ...
                'InitialGuess', [], ...
                'OptimalSequence', [], ...
                'PredictedState', [], ...
                'ExitFlag', NaN, ...
                'Cost', NaN, ...
                'SolverOutput', struct(), ...
                'GoalReached', false, ...
                'UsedFallback', false);
        end
    end
end
