classdef TestLiuSafeCorridorCartesian < matlab.unittest.TestCase
    %TESTLIUSAFECORRIDORCARTESIAN Stage 1 Cartesian generator tests.

    methods (TestClassSetup)
        function addClassFolder(testCase)
            classFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture(classFolder));
        end
    end

    methods (Test)
        function emptyMapCachesOneCorridorPerSegment(testCase)
            map = binaryOccupancyMap(false(120, 120), 10);
            path = [2 2; 7 2; 9 7];
            generator = LiuSafeCorridor(0.3, 2);

            corridors = generator.generate(path, map);

            testCase.verifyEqual(numel(corridors), size(path, 1)-1);
            testCase.verifyEqual(generator.Path, path);
            testCase.verifyEqual(generator.CumulativeLength, ...
                [0; 5; 5+sqrt(29)], 'AbsTol', 1e-12);
            testCase.verifyEqual(generator.Corridors, corridors);
            for segmentIndex = 1:numel(corridors)
                corridor = corridors(segmentIndex);
                testCase.verifyEqual(corridor.SegmentIndex, segmentIndex);
                testCase.verifySize(corridor.Axy, [8 2]);
                testCase.verifySize(corridor.bxy, [8 1]);
                testCase.verifyEqual( ...
                    vecnorm(corridor.Axy, 2, 2), ones(8, 1), ...
                    'AbsTol', 1e-12);
                testCase.verifyLessThanOrEqual( ...
                    max(corridor.Axy*corridor.Endpoints ...
                        - corridor.bxy, [], 'all'), 1e-11);
                testCase.verifyGreaterThan( ...
                    polyarea(corridor.Vertices(1, :), ...
                        corridor.Vertices(2, :)), 0);
                [Axy, bxy] = generator.getCartesianConstraints(segmentIndex);
                testCase.verifyEqual(Axy, corridor.Axy);
                testCase.verifyEqual(bxy, corridor.bxy);
            end

            sharedWaypoint = path(2, :).';
            testCase.verifyLessThanOrEqual( ...
                max(corridors(1).Axy*sharedWaypoint-corridors(1).bxy), ...
                1e-11);
            testCase.verifyLessThanOrEqual( ...
                max(corridors(2).Axy*sharedWaypoint-corridors(2).bxy), ...
                1e-11);
        end

        function shiftedMapAndRotatedSegmentsPreserveCellClearance(testCase)
            map = binaryOccupancyMap(false(120, 160), 10);
            map.GridLocationInWorld = [-20 30];
            obstaclePoints = [-16 36; -12 34; -8 39; -5 35];
            setOccupancy(map, obstaclePoints, true(size(obstaclePoints, 1), 1));
            path = [-18 32; -13 32; -9 36; -5 38];
            generator = LiuSafeCorridor(0.25, 2.5);

            corridors = generator.generate(path, map);

            testCase.verifyEqual(numel(corridors), 3);
            testCase.verifyCellSeparation(map, corridors, generator.RobotRadius);
            for segmentIndex = 1:numel(corridors)
                testCase.verifyLessThanOrEqual( ...
                    max(corridors(segmentIndex).Axy ...
                        * corridors(segmentIndex).Endpoints ...
                        - corridors(segmentIndex).bxy, [], 'all'), 1e-10);
            end
        end

        function normalAdjustmentPreservesNearEndpointSeed(testCase)
            map = binaryOccupancyMap(12, 12, 10);
            obstaclePoints = [9.25 5.15; 9.45 5.25; 6.05 7.05];
            setOccupancy(map, obstaclePoints, true(size(obstaclePoints, 1), 1));
            path = [3 5; 9 5];
            generator = LiuSafeCorridor(0.2, 3);

            corridor = generator.generate(path, map);

            testCase.verifyGreaterThan(corridor.AdjustedPlaneCount, 0);
            testCase.verifyLessThanOrEqual( ...
                max(corridor.Axy*corridor.Endpoints-corridor.bxy, ...
                    [], 'all'), 1e-10);
            testCase.verifyCellSeparation(map, corridor, generator.RobotRadius);
        end

        function diagonalBoxCornerCellUsesEuclideanSelection(testCase)
            map = binaryOccupancyMap(false(500, 500), 100);
            obstaclePoint = [3.505 3.005];
            setOccupancy(map, obstaclePoint, true);
            path = [0.999 2.499; 2.999 2.499];
            generator = LiuSafeCorridor(0.1, 0.5);

            corridor = generator.generate(path, map);

            % The cell disk misses the raw box diagonally and therefore
            % does not create an obstacle plane. The four box and four map
            % rows still keep the final corridor at Euclidean clearance.
            testCase.verifySize(corridor.Axy, [8 2]);
            testCase.verifyCellSeparation(map, corridor, generator.RobotRadius);
        end

        function generationIsDeterministicAndDoesNotMutateInputs(testCase)
            map = binaryOccupancyMap(14, 16, 8);
            map.GridLocationInWorld = [-3 -4];
            obstaclePoints = [3 3; 5 7; 8 2];
            setOccupancy(map, obstaclePoints, true(size(obstaclePoints, 1), 1));
            path = [-1 0; 4 1; 10 6];
            pathBefore = path;
            mapBefore = testCase.mapSnapshot(map);
            generator = LiuSafeCorridor(0.2, 2);

            first = generator.generate(path, map);
            second = generator.generate(path, map);

            testCase.verifyEqual(second, first);
            testCase.verifyEqual(path, pathBefore);
            testCase.verifyEqual(testCase.mapSnapshot(map), mapBefore);
        end

        function largeWorldOffsetKeepsPositiveAreaCorridor(testCase)
            map = binaryOccupancyMap(false(500, 500), 100);
            map.GridLocationInWorld = [1e8 -2e8];
            path = [1e8+1 -2e8+2; 1e8+3 -2e8+2];
            generator = LiuSafeCorridor(0.1, 0.5);

            corridor = generator.generate(path, map);

            localVertices = corridor.Vertices-corridor.Vertices(:, 1);
            testCase.verifyGreaterThan( ...
                polyarea(localVertices(1, :), localVertices(2, :)), 0);
            testCase.verifyLessThanOrEqual( ...
                max(corridor.Axy*corridor.Endpoints-corridor.bxy, ...
                    [], 'all'), 1e-6);
        end

        function visualizationSmokeTest(testCase)
            oldVisibility = get(groot, 'DefaultFigureVisible');
            cleanup = onCleanup(@() set( ...
                groot, 'DefaultFigureVisible', oldVisibility));
            set(groot, 'DefaultFigureVisible', 'off');

            map = binaryOccupancyMap(10, 10, 10);
            generator = LiuSafeCorridor(0.2, 2);
            generator.generate([2 2; 8 2], map);
            ax = generator.plot();
            figureHandle = ancestor(ax, 'figure');
            figureCleanup = onCleanup(@() close(figureHandle));

            testCase.verifyTrue(isgraphics(ax, 'axes'));
            testCase.verifyEqual(ax.XLim, map.XWorldLimits, 'AbsTol', 1e-12);
            testCase.verifyEqual(ax.YLim, map.YWorldLimits, 'AbsTol', 1e-12);
        end
    end

    methods (Access = private)
        function verifyCellSeparation(testCase, map, corridors, robotRadius)
            [rows, columns] = find(occupancyMatrix(map));
            if isempty(rows)
                return
            end
            points = grid2world(map, [rows, columns]).';
            clearance = robotRadius + sqrt(2)/(2*map.Resolution);

            for corridor = corridors
                testCase.verifyEqual( ...
                    vecnorm(corridor.Axy, 2, 2), ...
                    ones(size(corridor.bxy)), 'AbsTol', 1e-12);
                vertices = corridor.Vertices;
                for pointIndex = 1:size(points, 2)
                    point = points(:, pointIndex);
                    [inside, onBoundary] = inpolygon( ...
                        point(1), point(2), vertices(1, :), vertices(2, :));
                    testCase.verifyFalse(inside || onBoundary);

                    minimumDistance = Inf;
                    for vertexIndex = 1:size(vertices, 2)
                        nextIndex = mod(vertexIndex, size(vertices, 2)) + 1;
                        edgeStart = vertices(:, vertexIndex);
                        edge = vertices(:, nextIndex) - edgeStart;
                        fraction = dot(point-edgeStart, edge)/dot(edge, edge);
                        fraction = min(1, max(0, fraction));
                        projection = edgeStart + fraction*edge;
                        minimumDistance = min( ...
                            minimumDistance, norm(point-projection));
                    end
                    testCase.verifyGreaterThanOrEqual( ...
                        minimumDistance, clearance-1e-10);
                end
            end
        end

        function snapshot = mapSnapshot(~, map)
            snapshot = struct( ...
                'Occupancy', occupancyMatrix(map), ...
                'Resolution', map.Resolution, ...
                'GridSize', map.GridSize, ...
                'GridLocationInWorld', map.GridLocationInWorld, ...
                'XWorldLimits', map.XWorldLimits, ...
                'YWorldLimits', map.YWorldLimits, ...
                'DefaultValue', map.DefaultValue);
        end
    end
end
