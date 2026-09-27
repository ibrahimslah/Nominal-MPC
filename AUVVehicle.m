classdef AUVVehicle < handle
    %AUVVEHICLE Nonlinear planar AUV simulation plant in Cartesian coordinates.
    %
    % State:
    %   x = [px; py; yaw; vx; vy; r]
    %
    % Input:
    %   u = [Fx; tau]
    %
    % Dynamics:
    %
    %   px_dot  = vx*cos(yaw) - vy*sin(yaw)
    %   py_dot  = vx*sin(yaw) + vy*cos(yaw)
    %   yaw_dot = r
    %
    %   vx_dot = (Fx + my*vy*r - dx*vx) / mx
    %
    %   vy_dot = (-mx*vx*r - dy*vy) / my
    %
    %   r_dot  = (tau + (mx-my)*vx*vy - dr*r) / Iz
    %
    % Discrete stochastic propagation:
    %
    %   xMean(k+1) = f_d(x(k),u(k))
    %
    %   w(k) ~ N(0,Q)
    %
    %   x(k+1) = xMean(k+1) + w(k)
    %
    % where Q = StateCovariance.
    %
    % State propagation is performed using ode45 and the Gaussian
    % process noise is added after the deterministic propagation.


    properties (SetAccess = private)

        % Actual stochastic plant state
        State

        % Sampling time
        SampleTime

        % Predicted state before adding process noise
        % E[x(k+1) | x(k),u(k)]
        StateMean

        % Covariance of additive Gaussian process noise
        %
        % w(k) ~ N(0,StateCovariance)
        StateCovariance

    end


    properties (SetAccess = private)

        mx
        my
        Iz

        dx
        dy
        dr

    end


    properties (Constant)

        nx = 6
        nu = 2

    end


    methods

        %% ================================================================
        % Constructor
        % =================================================================

        function obj = AUVVehicle(vehicleParameters)
            %AUVVEHICLE Construct nonlinear stochastic Cartesian AUV plant.
            %
            % vehicleParameters fields:
            %
            %   .mx
            %   .my
            %   .Iz
            %
            %   .dx
            %   .dy
            %   .dr
            %
            %   .SampleTime
            %   .InitialState
            %
            %   .StateCovariance
            %
            %       6x6 covariance matrix of the additive Gaussian
            %       process noise:
            %
            %           w(k) ~ N(0,StateCovariance)


            % -------------------------------------------------------------
            % Vehicle parameters
            % -------------------------------------------------------------

            obj.mx = vehicleParameters.mx;
            obj.my = vehicleParameters.my;
            obj.Iz = vehicleParameters.Iz;

            obj.dx = vehicleParameters.dx;
            obj.dy = vehicleParameters.dy;
            obj.dr = vehicleParameters.dr;


            % -------------------------------------------------------------
            % Simulation parameters
            % -------------------------------------------------------------

            obj.SampleTime = vehicleParameters.SampleTime;


            % -------------------------------------------------------------
            % Initial state
            % -------------------------------------------------------------

            obj.State = vehicleParameters.InitialState;

            % At initialization the state mean is the supplied initial state
            obj.StateMean = obj.State;


            % -------------------------------------------------------------
            % Gaussian process-noise covariance
            % -------------------------------------------------------------

            if isfield(vehicleParameters, 'StateCovariance')

                obj.StateCovariance = ...
                    vehicleParameters.StateCovariance;

            else

                obj.StateCovariance = zeros(obj.nx);

            end


            % Force numerical symmetry
            obj.StateCovariance = ...
                0.5 * ...
                (obj.StateCovariance + ...
                 obj.StateCovariance.');

        end


        %% ================================================================
        % Step
        % =================================================================

        function xNext = step(obj, input)
            %STEP Propagate the stochastic nonlinear AUV model.
            %
            % The model is
            %
            %   x(k+1) = f_d(x(k),u(k)) + w(k)
            %
            % where
            %
            %   w(k) ~ N(0,Q)
            %
            % and
            %
            %   Q = obj.StateCovariance
            %
            %
            % First:
            %
            %   StateMean = f_d(State,input)
            %
            % is calculated using ode45.
            %
            % Then:
            %
            %   State = StateMean + w
            %
            % where w is a new independent Gaussian realization at
            % every call to step().
            %
            % Therefore the process noise is temporally white.


            % -------------------------------------------------------------
            % Deterministic nonlinear prediction
            % -------------------------------------------------------------

            tspan = [0 obj.SampleTime];


            [~, xTrajectory] = ode45( ...
                @(t, x) obj.dynamics(t, x, input), ...
                tspan, ...
                obj.State);


            % Predicted value / conditional mean
            %
            % E[x(k+1) | x(k),u(k)]

            obj.StateMean = ...
                xTrajectory(end, :).';


            % -------------------------------------------------------------
            % Generate zero-mean Gaussian white process noise
            % -------------------------------------------------------------
            %
            % w(k) ~ N(0,Q)
            %
            % Using eig() avoids requiring the Statistics and
            % Machine Learning Toolbox.
            %
            % A new randn realization is generated at every step,
            % therefore w(k) is temporally independent / white.

            [V, D] = eig(obj.StateCovariance);

            eigenvalues = ...
                max(diag(D), 0);

            covarianceSquareRoot = ...
                V * diag(sqrt(eigenvalues));


            gaussianNoise = ...
                covarianceSquareRoot * ...
                randn(obj.nx, 1);


            % -------------------------------------------------------------
            % Actual stochastic next state
            % -------------------------------------------------------------
            %
            % x(k+1) = xMean(k+1) + w(k)

            obj.State = ...
                obj.StateMean + gaussianNoise;


            % Return actual stochastic state

            xNext = obj.State;

        end


        %% ================================================================
        % Continuous nonlinear dynamics
        % =================================================================

        function xdot = dynamics(obj, ~, x, input)
            %DYNAMICS Continuous-time nonlinear Cartesian AUV dynamics.


            % -------------------------------------------------------------
            % States
            % -------------------------------------------------------------

            yaw = x(3);

            vx = x(4);
            vy = x(5);
            r  = x(6);


            % -------------------------------------------------------------
            % Inputs
            % -------------------------------------------------------------
            %
            % Input is assumed to be valid.
            %
            % input = [Fx; tau]

            Fx  = input(1);
            tau = input(2);


            % -------------------------------------------------------------
            % Kinematics
            % -------------------------------------------------------------

            pxDot = ...
                vx*cos(yaw) ...
                - vy*sin(yaw);


            pyDot = ...
                vx*sin(yaw) ...
                + vy*cos(yaw);


            yawDot = r;


            % -------------------------------------------------------------
            % Dynamics
            % -------------------------------------------------------------

            vxDot = ...
                ( ...
                    Fx ...
                    + obj.my*vy*r ...
                    - obj.dx*vx ...
                ) ...
                / obj.mx;


            vyDot = ...
                ( ...
                    -obj.mx*vx*r ...
                    - obj.dy*vy ...
                ) ...
                / obj.my;


            rDot = ...
                ( ...
                    tau ...
                    + (obj.mx - obj.my)*vx*vy ...
                    - obj.dr*r ...
                ) ...
                / obj.Iz;


            % -------------------------------------------------------------
            % State derivative
            % -------------------------------------------------------------

            xdot = [
                pxDot
                pyDot
                yawDot
                vxDot
                vyDot
                rDot
            ];

        end


        %% ================================================================
        % Set state
        % =================================================================

        function setState(obj, state)
            %SETSTATE Set plant state directly.

            obj.State = state;

            obj.StateMean = state;

        end


        %% ================================================================
        % Reset
        % =================================================================

        function reset(obj, state)
            %RESET Reset plant to supplied state.

            obj.State = state;

            obj.StateMean = state;

        end

    end

end