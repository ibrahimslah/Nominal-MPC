classdef TestAUVPathModel < matlab.unittest.TestCase
    %TESTAUVPATHMODEL Regression coverage for path-frame conversion.

    methods (TestClassSetup)
        function addProjectRoot(testCase)
            rootFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture(rootFolder));
        end
    end

    methods (Test)
        function projectionUsesFiniteDistanceAndRawAlongTrack(testCase)
            model = testCase.makeModel();
            cartesian = [15; 3; 0; 1.2; -0.4; 0.3];

            [state, segment] = model.cartesianToPath(cartesian);

            testCase.verifyEqual(segment, 2);
            testCase.verifyEqual(state(1), 8.275280724, 'AbsTol', 1e-9);
            testCase.verifyEqual(state(3), -4.494385525, 'AbsTol', 1e-9);
            testCase.verifyEqual(state(4:6), cartesian(4:6));
        end

        function updateIsStableAndReportsSelectedSegment(testCase)
            model = testCase.makeModel();
            cartesian = [15; 3; 5*pi; 1.2; -0.4; 0.3];

            [firstState, firstSegment] = model.updateFromCartesian(cartesian);
            [secondState, secondSegment] = model.updateFromCartesian(cartesian);

            testCase.verifyEqual(secondSegment, firstSegment);
            testCase.verifyEqual(secondState, firstState, 'AbsTol', 1e-12);
            testCase.verifyEqual(model.getState(), secondState, 'AbsTol', 1e-12);
        end

        function sharedWaypointsUseOutgoingSegments(testCase)
            model = testCase.makeModel();
            path = testCase.path();
            cartesian1 = [path(2, :).'; 0; 0; 0; 0];
            cartesian2 = [path(3, :).'; 0; 0; 0; 0];

            priorState = [model.CumulativeLength(2) + 2; 0; 0; 0; 0; 0];
            model.updateFromCartesian(model.pathToCartesian(priorState, 2));

            [state1, segment1] = model.updateFromCartesian(cartesian1);
            [state2, segment2] = model.updateFromCartesian(cartesian2);

            testCase.verifyEqual(segment1, 2);
            testCase.verifyEqual(segment2, 3);
            testCase.verifyEqual(state1(1), 8.5, 'AbsTol', 1e-12);
            testCase.verifyEqual(state2(1), 8.5 + sqrt(100.25), 'AbsTol', 1e-12);
        end

        function lateralErrorHasBothSignsOnInteriorSegments(testCase)
            model = testCase.makeModel();
            for segment = 1:3
                s = model.CumulativeLength(segment) + ...
                    model.SegmentLength(segment)/2;
                [position, ~, ~, normal] = model.pathGeometry(s, segment);
                for lateralError = [-2.25 1.75]
                    cartesian = [position.' + lateralError*normal.'; 0; 0; 0; 0];
                    [state, selectedSegment] = ...
                        model.cartesianToPath(cartesian, segment);
                    testCase.verifyEqual(selectedSegment, segment);
                    testCase.verifyEqual(state(1), s, 'AbsTol', 1e-11);
                    testCase.verifyEqual(state(3), lateralError, 'AbsTol', 1e-11);
                end
            end
        end

        function endpointExtrapolationRoundTripsExactly(testCase)
            model = testCase.makeModel();
            endpointS = [-2, model.TotalLength + 3];
            frames = [1, 3];
            for k = 1:numel(endpointS)
                state = [endpointS(k); 0.4; -1.3; 2; -3; 4];
                cartesian = model.pathToCartesian(state, frames(k));
                [converted, segment] = ...
                    model.cartesianToPath(cartesian, frames(k));
                testCase.verifyEqual(segment, frames(k));
                testCase.verifyEqual(converted, state, 'AbsTol', 1e-11);
            end
        end

        function conversionPreservesVelocityAndWrapsHeading(testCase)
            model = testCase.makeModel();
            s = model.CumulativeLength(2) + 3;
            [position, heading, ~, normal] = model.pathGeometry(s, 2);
            cartesian = [position.' + 0.7*normal.'; heading + 7*pi; ...
                4.2; -1.1; 0.6];

            [state, segment] = model.cartesianToPath(cartesian, 2);

            testCase.verifyEqual(segment, 2);
            testCase.verifyEqual(state(2), pi, 'AbsTol', 1e-12);
            testCase.verifyEqual(state(3), 0.7, 'AbsTol', 1e-12);
            testCase.verifyEqual(state(4:6), cartesian(4:6));
        end

        function ordinaryGeometryAndInverseRemainDefaultCompatible(testCase)
            model = testCase.makeModel();
            s = model.CumulativeLength(3) + 4;
            state = [s; -0.3; 0.9; 1; 2; 3];

            [position, heading, tangent, normal, segment] = model.pathGeometry(s);
            cartesian = model.pathToCartesian(state);
            expectedPosition = position + state(3)*normal;

            testCase.verifyEqual(segment, 3);
            testCase.verifyEqual(tangent, model.SegmentTangent(3, :));
            testCase.verifyEqual(heading, model.SegmentHeading(3), 'AbsTol', 1e-12);
            testCase.verifyEqual(cartesian(1:2), expectedPosition.', 'AbsTol', 1e-12);
            testCase.verifyEqual(cartesian(3), ...
                atan2(sin(heading + state(2)), cos(heading + state(2))), ...
                'AbsTol', 1e-12);
        end

        function nearWaypointTransitionsBetweenFrames(testCase)
            model = testCase.makeModel();
            waypointS = model.CumulativeLength(2);
            delta = 1e-9;
            before = [waypointS-delta; 0; 0; 0; 0; 0];
            after = [waypointS+delta; 0; 0; 0; 0; 0];

            beforeCartesian = model.pathToCartesian(before, 1);
            afterCartesian = model.pathToCartesian(after, 2);
            [beforeState, beforeSegment] = model.cartesianToPath(beforeCartesian);
            [afterState, afterSegment] = model.cartesianToPath(afterCartesian);

            testCase.verifyEqual(beforeSegment, 1);
            testCase.verifyEqual(afterSegment, 2);
            testCase.verifyEqual(beforeState(1), before(1), 'AbsTol', 1e-11);
            testCase.verifyEqual(afterState(1), after(1), 'AbsTol', 1e-11);
        end

        function unrelatedCrossingTiesHonorHintThenFirstSegment(testCase)
            model = testCase.makeCrossingModel();
            crossing = [1; 1; 0; 0; 0; 0];

            [withoutHint, firstSegment] = model.cartesianToPath(crossing);
            [withHint, hintedSegment] = model.cartesianToPath(crossing, 3);

            testCase.verifyEqual(firstSegment, 1);
            testCase.verifyEqual(hintedSegment, 3);
            testCase.verifyEqual(withoutHint(3), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(withHint(3), 0, 'AbsTol', 1e-12);
        end

        function segmentHintSearchRemainsLocal(testCase)
            model = testCase.makeModel();
            cartesian = [20; 9.81818181818182; 0; 0; 0; 0];

            [~, segment] = model.cartesianToPath(cartesian, 1);

            testCase.verifyEqual(segment, 2);
        end

        function explicitFramesRejectInvalidSegments(testCase)
            model = testCase.makeModel();
            state = [1; 0; 0; 0; 0; 0];

            identifier = 'AUVPathModel:InvalidFrameSegment';
            testCase.verifyError(@() model.pathToCartesian(state, 0), identifier);
            testCase.verifyError(@() model.pathToCartesian(state, 1.5), identifier);
            testCase.verifyError(@() model.pathGeometry(1, 4), identifier);
            testCase.verifyError(@() model.pathGeometry(1, []), identifier);
            testCase.verifyError(@() model.pathGeometry(1, 1+1i), identifier);
        end
    end

    methods (Static, Access = private)
        function model = makeModel()
            parameters = struct('SampleTime', 0.1);
            model = AUVPathModel(TestAUVPathModel.path(), parameters);
        end

        function model = makeCrossingModel()
            parameters = struct('SampleTime', 0.1);
            path = [0 0; 2 2; 0 2; 2 0];
            model = AUVPathModel(path, parameters);
        end

        function path = path()
            path = [2 3; 10.5 3; 10 13; 32 6];
        end

    end
end
