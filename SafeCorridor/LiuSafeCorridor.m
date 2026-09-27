classdef LiuSafeCorridor < handle
    %LIUSAFECORRIDOR Planar safe corridors around a piecewise-linear path.
    %
    % Implements the planar specialization of the ellipsoid-and-halfspace
    % construction in Section II-B of Liu et al., "Planning Dynamically
    % Feasible Trajectories for Quadrotors using Safe Flight Corridors in
    % 3-D Complex Environments". Occupied grid cells are conservatively
    % represented by circumscribed disks with radius
    % sqrt(2)/(2*map.Resolution).
    %
    % Public workflow:
    %   generator = LiuSafeCorridor(robotRadius, boundingBoxMargin);
    %   corridors = generator.generate(path, map);
    %   [Axy, bxy] = generator.getCartesianConstraints(segmentIndex);
    %   [Apath, bpath] = generator.getPathConstraints(segmentIndex);
    %   [Apath, bpath, segmentIndex] = ...
    %       generator.getPathConstraintsAtS(s);
    %   ax = generator.plot();
    %
    % robotRadius and boundingBoxMargin are in metres. path is an N-by-2
    % matrix of [x y] waypoints in metres, and map is a
    % binaryOccupancyMap. The caller supplies finite, nonzero-length path
    % segments with sufficient obstacle and map-boundary clearance, and a
    % boundingBoxMargin greater than robotRadius plus the occupied-cell
    % disk radius. Inputs are assumed correct and are not validated.
    %
    % Cartesian constraints use Axy*[x;y] <= bxy. Path constraints use the
    % state [s;ePsi;eY;vx;vy;r], where s is global cumulative progress and
    % positive eY points along the segment's left normal. For segment i,
    %
    %   [x;y] = p_i + t_i*(s-S_i) + n_i*eY.
    %
    % Apath has six columns, with nonzero entries only in columns 1 and 3,
    % and contains only converted corridor faces; SLimits are metadata and
    % are not appended as inequalities.
    %
    % Each generated corridor contains:
    %   SegmentIndex       path-segment index
    %   Endpoints          2-by-2 Cartesian endpoint columns
    %   SLimits            1-by-2 cumulative-progress interval
    %   Axy, bxy           m-by-2 and m-by-1 Cartesian inequalities
    %   Apath, bpath       m-by-6 and m-by-1 path inequalities
    %   Vertices           2-by-V Cartesian polygon vertices
    %   Ellipse            center, axes, semi-axes, E, and world-to-unit map
    %   AdjustedPlaneCount number of obstacle planes rotated for clearance
    %
    % getPathConstraintsAtS starts from half-open segment intervals, using
    % the outgoing segment at an internal knot and the first or final
    % segment outside the path interval. It selects the immediate next
    % corridor early when the current centerline point is inside it. Its
    % faces are returned in the path frame of the segment containing s.
    % Call generate before any getter or plot. The input path and map are
    % not modified. Navigation Toolbox is required for binaryOccupancyMap;
    % Optimization Toolbox is not used.

    properties (SetAccess = private)
        %ROBOTRADIUS Circular vehicle radius in metres.
        RobotRadius

        %BOUNDINGBOXMARGIN Liu local-box margin in metres.
        BoundingBoxMargin

        %PATH Last successfully generated N-by-2 waypoint path.
        Path = zeros(0, 2)

        %CUMULATIVELENGTH Arc length at every waypoint.
        CumulativeLength = zeros(0, 1)

        %CORRIDORS Cached Cartesian and path-coordinate corridor structs.
        Corridors = struct([])
    end

    properties (Access = private)
        %MAPSNAPSHOT Independent map copy used only for visualization.
        MapSnapshot = []
    end

    methods
        function obj = LiuSafeCorridor(robotRadius, boundingBoxMargin)
            %LIUSAFECORRIDOR Configure radius and local-box margin in metres.
            obj.RobotRadius = double(robotRadius);
            obj.BoundingBoxMargin = double(boundingBoxMargin);
        end

        function corridors = generate(obj, path, map)
            %GENERATE Build and cache one corridor per path segment.
            %
            % corridors = generate(path,map) accepts an N-by-2 waypoint path
            % and binaryOccupancyMap and returns a 1-by-(N-1) struct array.
            % The cache is replaced atomically: if any segment fails, the
            % most recent successful Path, CumulativeLength, Corridors, and
            % visualization map remain unchanged.
            [occupiedPoints, mapA, mapB, cellRadius] = obj.prepareMapGeometry(map);
            clearance = obj.RobotRadius + cellRadius;

            segmentVectors = diff(path, 1, 1);
            segmentLengths = vecnorm(segmentVectors, 2, 2);
            cumulativeLength = [0; cumsum(segmentLengths)];
            numberOfSegments = size(path, 1) - 1;

            emptyEllipse = struct( ...
                'Center', zeros(2, 1), ...
                'Axes', eye(2), ...
                'SemiAxes', zeros(2, 1), ...
                'E', zeros(2), ...
                'WorldToUnit', zeros(2));
            template = struct( ...
                'SegmentIndex', [], ...
                'Endpoints', zeros(2, 2), ...
                'SLimits', zeros(1, 2), ...
                'Axy', zeros(0, 2), ...
                'bxy', zeros(0, 1), ...
                'Apath', zeros(0, 6), ...
                'bpath', zeros(0, 1), ...
                'Vertices', zeros(2, 0), ...
                'Ellipse', emptyEllipse, ...
                'AdjustedPlaneCount', 0);
            newCorridors = repmat(template, 1, numberOfSegments);

            for segmentIndex = 1:numberOfSegments
                endpoints = path(segmentIndex:segmentIndex+1, :).';
                [boxA, boxB, box] = obj.buildSegmentBox(endpoints);
                localObstacles = obj.selectLocalObstacles( ...
                    occupiedPoints, box, cellRadius);

                ellipse = obj.fitEllipse(endpoints, localObstacles);
                [obstacleA, obstacleB, adjustedPlaneCount] = ...
                    obj.buildObstacleHalfspaces( ...
                        endpoints, ellipse, localObstacles, clearance);

                % Obstacle rows already include robot and cell clearance.
                % Shrinking the artificial box by RobotRadius guarantees
                % clearance from cells omitted by the local search, whose
                % circumscribed disks were included in the selection test.
                Axy = [obstacleA; boxA; mapA];
                bxy = [obstacleB; ...
                    boxB - obj.RobotRadius; ...
                    mapB - obj.RobotRadius];

                tangent = segmentVectors(segmentIndex, :).'/ ...
                    segmentLengths(segmentIndex);
                normal = [-tangent(2); tangent(1)];
                [Apath, bpath] = obj.cartesianToPathConstraints( ...
                    Axy, bxy, endpoints(:, 1), tangent, normal, ...
                    cumulativeLength(segmentIndex));

                vertices = obj.halfspacesToVertices( ...
                    Axy, bxy, map.XWorldLimits, map.YWorldLimits);

                newCorridors(segmentIndex).SegmentIndex = segmentIndex;
                newCorridors(segmentIndex).Endpoints = endpoints;
                newCorridors(segmentIndex).SLimits = ...
                    cumulativeLength(segmentIndex:segmentIndex+1).';
                newCorridors(segmentIndex).Axy = Axy;
                newCorridors(segmentIndex).bxy = bxy;
                newCorridors(segmentIndex).Apath = Apath;
                newCorridors(segmentIndex).bpath = bpath;
                newCorridors(segmentIndex).Vertices = vertices;
                newCorridors(segmentIndex).Ellipse = ellipse;
                newCorridors(segmentIndex).AdjustedPlaneCount = ...
                    adjustedPlaneCount;
            end

            % Commit only after every segment and the independent map copy
            % have been constructed successfully.
            mapSnapshot = copy(map);
            obj.Path = path;
            obj.CumulativeLength = cumulativeLength;
            obj.Corridors = newCorridors;
            obj.MapSnapshot = mapSnapshot;
            corridors = newCorridors;
        end

        function [Axy, bxy] = getCartesianConstraints(obj, segmentIndex)
            %GETCARTESIANCONSTRAINTS Return Axy*[x;y] <= bxy for a segment.
            Axy = obj.Corridors(segmentIndex).Axy;
            bxy = obj.Corridors(segmentIndex).bxy;
        end

        function [Apath, bpath] = getPathConstraints(obj, segmentIndex)
            %GETPATHCONSTRAINTS Return Apath*xPath <= bpath for one segment.
            Apath = obj.Corridors(segmentIndex).Apath;
            bpath = obj.Corridors(segmentIndex).bpath;
        end

        function [Apath, bpath, segmentIndex] = getPathConstraintsAtS(obj, s)
            %GETPATHCONSTRAINTSATS Select path constraints by cumulative s.
            % Uses the centerline point at s to select the next corridor
            % early when all its Cartesian halfspaces contain that point.
            % When selecting early, expresses the next corridor's faces in
            % the current segment's path frame. segmentIndex still names
            % the selected corridor.
            segmentIndex = find( s < (obj.CumulativeLength(2:end)), 1, 'first');
            if isempty(segmentIndex)
                segmentIndex = numel(obj.Corridors);
            end

            if segmentIndex < numel(obj.Corridors)
                segmentStart = obj.Path(segmentIndex, :).';
                segmentVector = obj.Path(segmentIndex+1, :).'-segmentStart;
                tangent = segmentVector/norm(segmentVector);
                position = segmentStart ...
                    + tangent*(s-obj.CumulativeLength(segmentIndex));

                nextCorridor = obj.Corridors(segmentIndex+1);
                tolerance = obj.scaledTolerance([position; nextCorridor.bxy]);
                insideNextCorridor = true;
                for faceIndex = 1:size(nextCorridor.Axy, 1)
                    if nextCorridor.Axy(faceIndex, :)*position ...
                            > nextCorridor.bxy(faceIndex)+tolerance
                        insideNextCorridor = false;
                        break
                    end
                end
                if insideNextCorridor
                    normal = [-tangent(2); tangent(1)];
                    [Apath, bpath] = obj.cartesianToPathConstraints( ...
                        nextCorridor.Axy, nextCorridor.bxy, ...
                        segmentStart, tangent, normal, ...
                        obj.CumulativeLength(segmentIndex));
                    segmentIndex = segmentIndex+1;
                    return
                end
            end

            Apath = obj.Corridors(segmentIndex).Apath;
            bpath = obj.Corridors(segmentIndex).bpath;
        end

        function ax = plot(obj)
            %PLOT Create a figure showing the cached map, path, and corridors.
            figureHandle = figure('Name', 'Liu 2-D safe corridors');
            ax = axes('Parent', figureHandle);
            show(obj.MapSnapshot, 'Parent', ax);
            hold(ax, 'on');

            colors = lines(max(1, numel(obj.Corridors)));
            for segmentIndex = 1:numel(obj.Corridors)
                vertices = obj.Corridors(segmentIndex).Vertices;
                patch(ax, vertices(1, :), vertices(2, :), ...
                    colors(segmentIndex, :), ...
                    'FaceAlpha', 0.18, ...
                    'EdgeColor', colors(segmentIndex, :), ...
                    'LineWidth', 1.2, ...
                    'DisplayName', sprintf('Safe corridor %d', segmentIndex));
            end
            plot(ax, obj.Path(:, 1), obj.Path(:, 2), 'r.-', ...
                'LineWidth', 1.6, 'MarkerSize', 15, ...
                'DisplayName', 'Reference path');

            axis(ax, 'equal');
            xlim(ax, obj.MapSnapshot.XWorldLimits);
            ylim(ax, obj.MapSnapshot.YWorldLimits);
            xlabel(ax, 'x (m)');
            ylabel(ax, 'y (m)');
            title(ax, 'Liu safe corridors and reference path');
            hold(ax, 'off');
        end
    end

    methods (Access = private)
        function [points, A, b, cellRadius] = prepareMapGeometry(~, map)
            occupied = occupancyMatrix(map);
            [rows, columns] = find(occupied);
            if isempty(rows)
                points = zeros(2, 0);
            else
                points = grid2world(map, [rows, columns]).';
            end

            A = [1 0; -1 0; 0 1; 0 -1];
            b = [map.XWorldLimits(2); -map.XWorldLimits(1); ...
                 map.YWorldLimits(2); -map.YWorldLimits(1)];
            cellRadius = sqrt(2)/(2*map.Resolution);
        end

        function [A, b, box] = buildSegmentBox(obj, endpoints)
            displacement = endpoints(:, 2) - endpoints(:, 1);
            segmentLength = norm(displacement);
            tangent = displacement/segmentLength;
            normal = [-tangent(2); tangent(1)];
            center = mean(endpoints, 2);

            A = [tangent.'; -tangent.'; normal.'; -normal.'];
            halfLengths = [segmentLength/2 + obj.BoundingBoxMargin; ...
                           segmentLength/2 + obj.BoundingBoxMargin; ...
                           obj.BoundingBoxMargin; ...
                           obj.BoundingBoxMargin];
            b = A*center + halfLengths;
            box = struct( ...
                'Center', center, ...
                'Axes', [tangent, normal], ...
                'HalfLengths', [halfLengths(1); halfLengths(3)]);
        end

        function localPoints = selectLocalObstacles(obj, points, box, cellRadius)
            if isempty(points)
                localPoints = zeros(2, 0);
                return
            end
            localPointsInBoxFrame = box.Axes.'*(points - box.Center);
            tolerance = obj.scaledTolerance( ...
                [box.HalfLengths; localPointsInBoxFrame(:); cellRadius]);
            excess = max( ...
                abs(localPointsInBoxFrame) - box.HalfLengths, 0);
            intersectsBox = vecnorm(excess, 2, 1) ...
                <= cellRadius + tolerance;
            localPoints = points(:, intersectsBox);
        end

        function ellipse = fitEllipse(~, endpoints, points)
            displacement = endpoints(:, 2) - endpoints(:, 1);
            segmentLength = norm(displacement);
            majorAxis = segmentLength/2;
            tangent = displacement/segmentLength;
            normal = [-tangent(2); tangent(1)];
            frame = [tangent, normal];
            center = mean(endpoints, 2);

            localPoints = frame.'*(points - center);
            active = abs(localPoints(1, :)) < majorAxis;
            normalizedAlong = localPoints(1, active)/majorAxis;
            denominators = sqrt( ...
                (1 - normalizedAlong).*(1 + normalizedAlong));
            minorBounds = abs(localPoints(2, active))./denominators;
            minorAxis = min([majorAxis, minorBounds]);

            semiAxes = [majorAxis; minorAxis];
            ellipse = struct( ...
                'Center', center, ...
                'Axes', frame, ...
                'SemiAxes', semiAxes, ...
                'E', frame*diag(semiAxes)*frame.', ...
                'WorldToUnit', diag(1./semiAxes)*frame.');
        end

        function [A, b, adjusted] = buildObstacleHalfspaces( ...
                obj, endpoints, ellipse, points, clearance)
            numberOfPoints = size(points, 2);
            A = zeros(numberOfPoints, 2);
            b = zeros(numberOfPoints, 1);
            adjusted = 0;
            if numberOfPoints == 0
                return
            end

            F = ellipse.WorldToUnit;
            unitPoints = F*(points - ellipse.Center);
            metric = sum(unitPoints.^2, 1);
            remaining = 1:numberOfPoints;
            planeCount = 0;
            relativePoints = points - ellipse.Center;
            relativeEndpoints = endpoints - ellipse.Center;
            tolerance = obj.scaledTolerance( ...
                [relativePoints(:); relativeEndpoints(:); ...
                 ellipse.SemiAxes; clearance]);

            while ~isempty(remaining)
                [~, localIndex] = min(metric(remaining));
                contactIndex = remaining(localIndex);
                contact = points(:, contactIndex);

                gradient = F.'*unitPoints(:, contactIndex);
                normal = gradient/norm(gradient);

                endpointGaps = normal.'*(contact - endpoints);
                if min(endpointGaps) < clearance - tolerance
                    normal = obj.adjustNormal( ...
                        normal, contact, endpoints, clearance, tolerance);
                    adjusted = adjusted + 1;
                end

                planeCount = planeCount + 1;
                A(planeCount, :) = normal.';
                b(planeCount) = normal.'*contact - clearance;

                % Recompute pruning with an adjusted normal. At least the
                % selected contact is removed, guaranteeing termination.
                remove = normal.'*(points(:, remaining) - contact) >= 0;
                remove(localIndex) = true;
                remaining = remaining(~remove);
            end

            A = A(1:planeCount, :);
            b = b(1:planeCount);
        end

        function normal = adjustNormal(~, original, contact, endpoints, ...
                clearance, tolerance)
            endpointVectors = contact - endpoints;
            candidates = zeros(2, 0);

            for endpointIndex = 1:2
                vector = endpointVectors(:, endpointIndex);
                distance = norm(vector);
                ratio = min(1, max(-1, clearance/distance));
                centerAngle = atan2(vector(2), vector(1));
                spread = acos(ratio);
                angles = centerAngle + [-spread, spread];
                candidates = [candidates, ...
                    [cos(angles); sin(angles)]]; %#ok<AGROW>
            end

            feasible = all( ...
                candidates.'*endpointVectors >= clearance - tolerance, 2);
            candidates = candidates(:, feasible);

            [~, choice] = max(original.'*candidates);
            normal = candidates(:, choice);
        end

        function vertices = halfspacesToVertices( obj, A, b, xLimits, yLimits)
            origin = [mean(xLimits); mean(yLimits)];
            localXLimits = xLimits - origin(1);
            localYLimits = yLimits - origin(2);
            localB = b - A*origin;
            vertices = [localXLimits([1 2 2 1]); ...
                        localYLimits([1 1 2 2])];
            tolerance = obj.scaledTolerance( ...
                [localB; localXLimits(:); localYLimits(:)]);

            for face = 1:size(A, 1)
                if isempty(vertices)
                    return
                end

                output = zeros(2, 0);
                startPoint = vertices(:, end);
                startResidual = A(face, :)*startPoint - localB(face);
                startInside = startResidual <= tolerance;

                for vertexIndex = 1:size(vertices, 2)
                    endPoint = vertices(:, vertexIndex);
                    endResidual = A(face, :)*endPoint - localB(face);
                    endInside = endResidual <= tolerance;

                    if startInside ~= endInside
                        fraction = startResidual ...
                            / (startResidual - endResidual);
                        output(:, end+1) = startPoint ...
                            + fraction*(endPoint - startPoint); %#ok<AGROW>
                    end
                    if endInside
                        output(:, end+1) = endPoint; %#ok<AGROW>
                    end

                    startPoint = endPoint;
                    startResidual = endResidual;
                    startInside = endInside;
                end
                vertices = output;
            end
            vertices = vertices + origin;
        end

        function [Apath, bpath] = cartesianToPathConstraints( ~, Axy, bxy, segmentStart, tangent, normal, segmentStartS)
            numberOfConstraints = size(Axy, 1);
            Apath = zeros(numberOfConstraints, 6);
            Apath(:, 1) = Axy*tangent;
            Apath(:, 3) = Axy*normal;
            bpath = bxy - Axy*segmentStart + Apath(:, 1)*segmentStartS;
        end

        function tolerance = scaledTolerance(~, values)
            scale = max([1; abs(values(:))]);
            tolerance = 128*eps(scale);
        end
    end
end
