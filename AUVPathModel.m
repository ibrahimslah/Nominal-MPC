classdef AUVPathModel < handle
    %AUVPATHMODEL Stateful path-coordinate model for planar AUV MPC.
    %
    % Path state:
    %
    %   x = [s;
    %        ePsi;
    %        eY;
    %        vx;
    %        vy;
    %        r]
    %
    % Input:
    %
    %   u = [Fx;
    %        tau]
    %
    % Main use:
    %
    %   updateFromCartesian()   Synchronize with plant state
    %   getState()              Get current path state
    %   step()                  Propagate internal path state
    %   nonlinearDiscrete()     Nonlinear prediction model
    %   linearizedDiscrete()    Linearized MPC model
    %   cartesianToPath()       Coordinate transformation


    properties (SetAccess = private)

        State = zeros(6,1)

        Path

        Parameters

        SampleTime

        SegmentVector
        SegmentLength
        SegmentTangent
        SegmentNormal
        SegmentHeading

        CumulativeLength
        TotalLength

        ActiveSegment = 1

    end


    properties (Constant)

        nx = 6
        nu = 2

    end


    methods

        % ================================================================
        %% Constructor
        % =================================================================

        function obj = AUVPathModel(path, vehicleParameters, initialState)

            obj.Path = path;

            obj.Parameters = vehicleParameters;

            obj.SampleTime = ...
                vehicleParameters.SampleTime;


            % -------------------------------------------------------------
            % Path geometry
            % -------------------------------------------------------------

            obj.SegmentVector = ...
                diff(path,1,1);


            obj.SegmentLength = ...
                sqrt(sum(obj.SegmentVector.^2,2));


            obj.SegmentTangent = ...
                obj.SegmentVector ./ obj.SegmentLength;


            obj.SegmentNormal = [
                -obj.SegmentTangent(:,2), ...
                 obj.SegmentTangent(:,1)
            ];


            obj.SegmentHeading = atan2( ...
                obj.SegmentTangent(:,2), ...
                obj.SegmentTangent(:,1));


            obj.CumulativeLength = [
                0
                cumsum(obj.SegmentLength)
            ];


            obj.TotalLength = ...
                obj.CumulativeLength(end);


            % -------------------------------------------------------------
            % Initial path-coordinate state
            % -------------------------------------------------------------

            if nargin >= 3

                obj.State = initialState;

            else

                obj.State = zeros(6,1);

            end

        end


        % ================================================================
        %% Get internal state
        % =================================================================

        function x = getState(obj)

            x = obj.State;

        end


        % ================================================================
        %% Set internal path state
        % =================================================================

        function setState(obj, x)

            obj.State = x;

        end


        % ================================================================
        %% Synchronize path state with Cartesian plant
        % =================================================================

        function xPath = updateFromCartesian(obj, xCartesian)
            %UPDATEFROMCARTESIAN
            %
            % Converts the current Cartesian plant state into path
            % coordinates and stores it as the internal State.
            %
            % Cartesian state:
            %
            %   [px;
            %    py;
            %    psi;
            %    vx;
            %    vy;
            %    r]
            %
            % Path state:
            %
            %   [s;
            %    ePsi;
            %    eY;
            %    vx;
            %    vy;
            %    r]


            [xPath, segmentIndex] = ...
                obj.cartesianToPath( ...
                    xCartesian, ...
                    obj.ActiveSegment);


            obj.State = xPath;

            obj.ActiveSegment = segmentIndex;

        end


        %% ================================================================
        % Advance internal state
        % =================================================================

        function xNext = step(obj, u)
            %STEP Propagate internal path-coordinate state by one sample.

            obj.State = ...
                obj.nonlinearDiscrete( ...
                    obj.State, ...
                    u);

            xNext = obj.State;

        end


        %% ================================================================
        % Nonlinear discrete model
        % =================================================================

        function xNext = nonlinearDiscrete(obj, x, u)

            s    = x(1);
            ePsi = x(2);
            eY   = x(3);
            vx   = x(4);
            vy   = x(5);
            r    = x(6);


            Fx  = u(1);
            tau = u(2);


            p = obj.Parameters;

            mx = p.mx;
            my = p.my;
            Iz = p.Iz;

            dx = p.dx;
            dy = p.dy;
            dr = p.dr;

            Ts = obj.SampleTime;


            % -------------------------------------------------------------
            % Piecewise-linear path:
            %
            % kappa = 0
            % -------------------------------------------------------------

            C = cos(ePsi);
            S = sin(ePsi);


            Vt = ...
                vx*C - vy*S;


            Vn = ...
                vx*S + vy*C;


            % -------------------------------------------------------------
            % Path-coordinate state update
            % -------------------------------------------------------------

            sNext = ...
                s + Ts*Vt;


            ePsiNext = ...
                ePsi + Ts*r;


            eYNext = ...
                eY + Ts*Vn;


            % -------------------------------------------------------------
            % Dynamic state update
            % -------------------------------------------------------------

            vxNext = ...
                vx ...
                + Ts/mx * ...
                ( ...
                    Fx ...
                    + my*vy*r ...
                    - dx*vx ...
                );


            vyNext = ...
                vy ...
                + Ts/my * ...
                ( ...
                    -mx*vx*r ...
                    - dy*vy ...
                );


            rNext = ...
                r ...
                + Ts/Iz * ...
                ( ...
                    tau ...
                    + (mx-my)*vx*vy ...
                    - dr*r ...
                );


            xNext = [
                sNext
                ePsiNext
                eYNext
                vxNext
                vyNext
                rNext
            ];

        end


        %% ================================================================
        % Discrete linearization
        % =================================================================

        function [Ad, Bd, cd] = linearizedDiscrete(obj, x, u)

            ePsi = x(2);

            vx = x(4);
            vy = x(5);
            r  = x(6);


            p = obj.Parameters;

            mx = p.mx;
            my = p.my;
            Iz = p.Iz;

            dx = p.dx;
            dy = p.dy;
            dr = p.dr;

            Ts = obj.SampleTime;


            C = cos(ePsi);
            S = sin(ePsi);


            Vt = ...
                vx*C - vy*S;

            Vn = ...
                vx*S + vy*C;


            % -------------------------------------------------------------
            % State Jacobian
            % -------------------------------------------------------------

            Ad = zeros(6,6);


            % s+
            Ad(1,1) = 1;

            Ad(1,2) = ...
                -Ts*Vn;

            Ad(1,4) = ...
                Ts*C;

            Ad(1,5) = ...
                -Ts*S;


            % ePsi+
            Ad(2,2) = 1;

            Ad(2,6) = Ts;


            % eY+
            Ad(3,2) = ...
                Ts*Vt;

            Ad(3,3) = 1;

            Ad(3,4) = ...
                Ts*S;

            Ad(3,5) = ...
                Ts*C;


            % vx+
            Ad(4,4) = ...
                1 - Ts*dx/mx;

            Ad(4,5) = ...
                Ts*my*r/mx;

            Ad(4,6) = ...
                Ts*my*vy/mx;


            % vy+
            Ad(5,4) = ...
                -Ts*mx*r/my;

            Ad(5,5) = ...
                1 - Ts*dy/my;

            Ad(5,6) = ...
                -Ts*mx*vx/my;


            % r+
            Ad(6,4) = ...
                Ts*(mx-my)*vy/Iz;

            Ad(6,5) = ...
                Ts*(mx-my)*vx/Iz;

            Ad(6,6) = ...
                1 - Ts*dr/Iz;


            % -------------------------------------------------------------
            % Input Jacobian
            % -------------------------------------------------------------

            Bd = zeros(6,2);

            Bd(4,1) = ...
                Ts/mx;

            Bd(6,2) = ...
                Ts/Iz;


            % -------------------------------------------------------------
            % Affine term
            % -------------------------------------------------------------

            xNext = ...
                obj.nonlinearDiscrete(x,u);


            cd = ...
                xNext ...
                - Ad*x ...
                - Bd*u;

        end


        %% ================================================================
        % Cartesian -> path coordinates
        % =================================================================

        function [xPath, segmentIndex] = cartesianToPath(obj, xCartesian, segmentHint)


            px  = xCartesian(1);
            py  = xCartesian(2);
            psi = xCartesian(3);

            position = [px py];


            numberOfSegments = ...
                length(obj.SegmentLength);


            % -------------------------------------------------------------
            % Search region
            % -------------------------------------------------------------

            if nargin < 3

                candidateSegments = ...
                    1:numberOfSegments;

            else

                firstSegment = ...
                    max(1,segmentHint-1);

                lastSegment = ...
                    min( ...
                        numberOfSegments, ...
                        segmentHint+1);

                candidateSegments = ...
                    firstSegment:lastSegment;

            end


            minimumDistanceSquared = Inf;

            segmentIndex = ...
                candidateSegments(1);

            bestAlongSegment = 0;

            bestProjection = ...
                obj.Path(segmentIndex,:);


            % -------------------------------------------------------------
            % Closest projection
            % -------------------------------------------------------------

            for i = candidateSegments

                segmentStart = ...
                    obj.Path(i,:);

                tangent = ...
                    obj.SegmentTangent(i,:);

                L = ...
                    obj.SegmentLength(i);


                localS = ...
                    dot( ...
                        position-segmentStart, ...
                        tangent);


                localS = ...
                    min(max(localS,0),L);


                projection = ...
                    segmentStart ...
                    + localS*tangent;


                errorVector = ...
                    position-projection;


                distanceSquared = ...
                    dot(errorVector,errorVector);


                if distanceSquared < ...
                        minimumDistanceSquared

                    minimumDistanceSquared = ...
                        distanceSquared;

                    segmentIndex = i;

                    bestAlongSegment = ...
                        localS;

                    bestProjection = ...
                        projection;

                end

            end


            % -------------------------------------------------------------
            % Global path distance
            % -------------------------------------------------------------

            s = ...
                obj.CumulativeLength(segmentIndex) ...
                + bestAlongSegment;


            % -------------------------------------------------------------
            % Signed lateral error
            % -------------------------------------------------------------

            normal = ...
                obj.SegmentNormal(segmentIndex,:);


            eY = ...
                dot( ...
                    position-bestProjection, ...
                    normal);


            % -------------------------------------------------------------
            % Heading error
            % -------------------------------------------------------------

            pathHeading = ...
                obj.SegmentHeading(segmentIndex);


            ePsi = ...
                atan2( ...
                    sin(psi-pathHeading), ...
                    cos(psi-pathHeading));


            % -------------------------------------------------------------
            % Result
            % -------------------------------------------------------------

            xPath = [
                s
                ePsi
                eY
                xCartesian(4)
                xCartesian(5)
                xCartesian(6)
            ];

        end


        % ================================================================
        %% Path -> Cartesian coordinates
        % =================================================================

        function xCartesian = pathToCartesian(obj, xPath)
            %PATHTOCARTESIAN
            %
            % Useful for:
            %
            %   - plotting predictions
            %   - debugging MPC
            %   - checking coordinate transformations


            s     = xPath(1);
            ePsi  = xPath(2);
            eY    = xPath(3);


            [referencePosition, pathHeading, ~, normal] = ...
                obj.pathGeometry(s);


            position = ...
                referencePosition ...
                + eY*normal;


            psi = ...
                obj.wrapAngle(pathHeading + ePsi);


            xCartesian = [
                position(1)
                position(2)
                psi
                xPath(4)
                xPath(5)
                xPath(6)
            ];

        end


        % ================================================================
        %% Get local path geometry from s
        % =================================================================

        function [position, heading, tangent, normal, segmentIndex] = pathGeometry(obj, s)
            %PATHGEOMETRY
            %
            % Return local geometry of the piecewise-linear path at
            % global arc-length coordinate s.
            %
            % The first and final segments are extrapolated if s lies
            % slightly outside [0, TotalLength].


            numberOfSegments = ...
                length(obj.SegmentLength);


            % -------------------------------------------------------------
            % Determine segment
            % -------------------------------------------------------------

            if s <= 0

                segmentIndex = 1;

                localS = s;


            elseif s >= obj.TotalLength

                segmentIndex = numberOfSegments;

                localS = obj.SegmentLength(end) + (s-obj.TotalLength);


            else

                segmentIndex = find( s < obj.CumulativeLength(2:end), 1, 'first');


                localS = s - obj.CumulativeLength(segmentIndex);

            end


            tangent = obj.SegmentTangent(segmentIndex,:);


            normal = obj.SegmentNormal(segmentIndex,:);


            heading = obj.SegmentHeading(segmentIndex);


            segmentStart = obj.Path(segmentIndex,:);


            position = segmentStart + localS*tangent;

        end


        % ================================================================
        %% Generate nominal path-following reference
        % =================================================================

        function xReference = referenceState(~, sReference, vxReference)
            %REFERENCESTATE
            %
            % Nominal Frenet-frame path-following state:
            %
            %   ePsi = 0
            %   eY   = 0
            %   vy   = 0
            %   r    = 0


            xReference = [
                sReference
                0
                0
                vxReference
                0
                0
            ];

        end

    end


    methods (Static, Access = private)

        % ================================================================
        %% Angle wrapping
        % =================================================================

        function angle = wrapAngle(angle)

            angle = ...
                atan2( ...
                    sin(angle), ...
                    cos(angle));

        end

    end

end