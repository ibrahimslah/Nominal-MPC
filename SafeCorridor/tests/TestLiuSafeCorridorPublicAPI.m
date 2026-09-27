classdef TestLiuSafeCorridorPublicAPI < matlab.unittest.TestCase
    %TESTLIUSAFECORRIDORPUBLICAPI Black-box tests for the public API.

    methods (TestClassSetup)
        function addClassFolder(testCase)
            classFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture(classFolder));
        end
    end

    methods (Test)
        function generatedCorridorsExposeCompletePublicContract(testCase)
            [generator, corridors, path] = testCase.generateShiftedFixture();
            expectedFields = { ...
                'SegmentIndex'; 'Endpoints'; 'SLimits'; 'Axy'; 'bxy'; ...
                'Apath'; 'bpath'; 'Vertices'; 'Ellipse'; ...
                'AdjustedPlaneCount'};
            expectedEllipseFields = { ...
                'Center'; 'Axes'; 'SemiAxes'; 'E'; 'WorldToUnit'};
            expectedCumulativeLength = [0; cumsum( ...
                vecnorm(diff(path, 1, 1), 2, 2))];

            testCase.verifyEqual(generator.RobotRadius, 0.25);
            testCase.verifyEqual(generator.BoundingBoxMargin, 2.5);
            testCase.verifyEqual(generator.Path, path);
            testCase.verifyEqual(generator.CumulativeLength, ...
                expectedCumulativeLength, 'AbsTol', 1e-12);
            testCase.verifyEqual(generator.Corridors, corridors);
            testCase.verifyTrue(isstruct(corridors));
            testCase.verifySize(corridors, [1, size(path, 1)-1]);
            testCase.verifyEqual(fieldnames(corridors), expectedFields);

            hasObstaclePlane = false;
            for segmentIndex = 1:numel(corridors)
                corridor = corridors(segmentIndex);
                constraintCount = size(corridor.Axy, 1);

                testCase.verifyEqual(corridor.SegmentIndex, segmentIndex);
                testCase.verifyEqual(corridor.Endpoints, ...
                    path(segmentIndex:segmentIndex+1, :).');
                testCase.verifyEqual(corridor.SLimits, ...
                    expectedCumulativeLength( ...
                    segmentIndex:segmentIndex+1).', 'AbsTol', 1e-12);
                testCase.verifyGreaterThanOrEqual(constraintCount, 8);
                testCase.verifySize(corridor.Axy, [constraintCount, 2]);
                testCase.verifySize(corridor.bxy, [constraintCount, 1]);
                testCase.verifySize(corridor.Apath, [constraintCount, 6]);
                testCase.verifySize(corridor.bpath, [constraintCount, 1]);
                constraintValues = [corridor.Axy(:); corridor.bxy(:); ...
                    corridor.Apath(:); corridor.bpath(:)];
                testCase.verifyTrue(all(isfinite(constraintValues)));
                testCase.verifyEqual( ...
                    vecnorm(corridor.Axy, 2, 2), ...
                    ones(constraintCount, 1), 'AbsTol', 1e-12);
                testCase.verifyLessThanOrEqual( ...
                    max(corridor.Axy*corridor.Endpoints ...
                    - corridor.bxy, [], 'all'), 1e-10);
                testCase.verifyEqual(corridor.Apath(:, [2 4 5 6]), ...
                    zeros(constraintCount, 4));
                testCase.verifyEqual(size(corridor.Vertices, 1), 2);
                testCase.verifyGreaterThanOrEqual( ...
                    size(corridor.Vertices, 2), 3);
                testCase.verifyTrue(all(isfinite(corridor.Vertices), 'all'));
                testCase.verifySize(corridor.AdjustedPlaneCount, [1 1]);
                testCase.verifyGreaterThanOrEqual( ...
                    corridor.AdjustedPlaneCount, 0);
                testCase.verifyEqual(corridor.AdjustedPlaneCount, ...
                    fix(corridor.AdjustedPlaneCount));

                [Axy, bxy] = generator.getCartesianConstraints(segmentIndex);
                [Apath, bpath] = generator.getPathConstraints(segmentIndex);
                testCase.verifyEqual(Axy, corridor.Axy);
                testCase.verifyEqual(bxy, corridor.bxy);
                testCase.verifyEqual(Apath, corridor.Apath);
                testCase.verifyEqual(bpath, corridor.bpath);

                ellipse = corridor.Ellipse;
                testCase.verifyEqual(fieldnames(ellipse), ...
                    expectedEllipseFields);
                testCase.verifySize(ellipse.Center, [2 1]);
                testCase.verifySize(ellipse.Axes, [2 2]);
                testCase.verifySize(ellipse.SemiAxes, [2 1]);
                testCase.verifySize(ellipse.E, [2 2]);
                testCase.verifySize(ellipse.WorldToUnit, [2 2]);
                ellipseValues = [ellipse.Center(:); ellipse.Axes(:); ...
                    ellipse.SemiAxes(:); ellipse.E(:); ...
                    ellipse.WorldToUnit(:)];
                testCase.verifyTrue(all(isfinite(ellipseValues)));
                testCase.verifyGreaterThan(ellipse.SemiAxes, zeros(2, 1));
                testCase.verifyEqual(ellipse.Center, ...
                    mean(corridor.Endpoints, 2), 'AbsTol', 1e-12);
                testCase.verifyEqual(ellipse.Axes.'*ellipse.Axes, eye(2), ...
                    'AbsTol', 1e-12);
                testCase.verifyEqual(ellipse.E, ...
                    ellipse.Axes*diag(ellipse.SemiAxes)*ellipse.Axes.', ...
                    'AbsTol', 1e-12);
                testCase.verifyEqual(ellipse.WorldToUnit ...
                    * ellipse.Axes*diag(ellipse.SemiAxes), eye(2), ...
                    'AbsTol', 1e-12);

                hasObstaclePlane = hasObstaclePlane || constraintCount > 8;
            end
            testCase.verifyTrue(hasObstaclePlane);
        end

        function pathResidualsMatchCartesianResiduals(testCase)
            [generator, corridors] = testCase.generateShiftedFixture();
            lateralOffsets = [-0.4, 0, 0.4];
            inactiveState = [0.23; 1.4; -0.6; 0.11];
            alternateInactiveState = [-0.81; -2.2; 0.7; -0.31];

            for segmentIndex = 1:numel(corridors)
                corridor = corridors(segmentIndex);
                [Axy, bxy] = ...
                    generator.getCartesianConstraints(segmentIndex);
                [Apath, bpath] = ...
                    generator.getPathConstraints(segmentIndex);
                displacement = corridor.Endpoints(:, 2) ...
                    - corridor.Endpoints(:, 1);
                tangent = displacement/norm(displacement);
                normal = [-tangent(2); tangent(1)];
                sSamples = [corridor.SLimits(1), ...
                    mean(corridor.SLimits), corridor.SLimits(2)];

                for s = sSamples
                    for eY = lateralOffsets
                        state = [s; inactiveState(1); eY; ...
                            inactiveState(2:4)];
                        position = corridor.Endpoints(:, 1) ...
                            + tangent*(s-corridor.SLimits(1)) ...
                            + normal*eY;

                        testCase.verifyEqual(Apath*state-bpath, ...
                            Axy*position-bxy, 'AbsTol', 1e-11);

                        alternateState = [s; alternateInactiveState(1); ...
                            eY; alternateInactiveState(2:4)];
                        testCase.verifyEqual( ...
                            Apath*alternateState-bpath, ...
                            Apath*state-bpath, 'AbsTol', 1e-12);
                    end
                end
            end
        end

        function onlineLookupUsesHalfOpenIntervalsAndClamps(testCase)
            [generator, corridors] = testCase.generateShiftedFixture();
            knots = generator.CumulativeLength;
            delta = min(diff(knots))*1e-8;
            queries = [ ...
                knots(1)-1, knots(1), ...
                knots(2)-delta, knots(2), knots(2)+delta, ...
                knots(3)-delta, knots(3), knots(3)+delta, ...
                knots(end)-delta, knots(end), knots(end)+1];
            expectedSegments = [1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3];

            testCase.verifyEqual(numel(corridors), 3);
            for queryIndex = 1:numel(queries)
                expectedSegment = expectedSegments(queryIndex);
                [expectedA, expectedB] = ...
                    generator.getPathConstraints(expectedSegment);
                [actualA, actualB, actualSegment] = ...
                    generator.getPathConstraintsAtS(queries(queryIndex));

                testCase.verifyEqual(actualSegment, expectedSegment);
                testCase.verifyEqual(actualA, expectedA);
                testCase.verifyEqual(actualB, expectedB);
            end
        end
    end

    methods (Access = private)
        function [generator, corridors, path] = generateShiftedFixture(~)
            map = binaryOccupancyMap(false(120, 160), 10);
            map.GridLocationInWorld = [-20 30];
            obstaclePoints = [-16 36; -12 34; -8 39; -5 35];
            setOccupancy(map, obstaclePoints, ...
                true(size(obstaclePoints, 1), 1));
            path = [-18 32; -13 32; -9 36; -5 38];
            generator = LiuSafeCorridor(0.25, 2.5);
            corridors = generator.generate(path, map);
        end
    end
end
