classdef TestPathFrameIntegration < matlab.unittest.TestCase
    %TESTPATHFRAMEINTEGRATION Preserve the measured frame through the MPC.

    methods (TestClassSetup)
        function addFolders(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'SafeCorridor')));
        end
    end

    methods (Test)
        function explicitFrameBelowKnotMatchesCartesianFaces(testCase)
            [model, corridors] = testCase.fixture();
            x = [15; 3; 0.3; 1; -0.2; 0.1];
            [z, frame] = model.updateFromCartesian(x);
            testCase.verifyEqual(frame, 2);
            testCase.verifyLessThan(z(1), model.CumulativeLength(2));
            [A, b, corridor] = corridors.getPathConstraintsAtS(z(1), frame);
            [Axy, bxy] = corridors.getCartesianConstraints(corridor);
            testCase.verifyEqual(A*z-b, Axy*x(1:2)-bxy, 'AbsTol', 1e-11);
            testCase.verifyEqual(model.pathToCartesian(z, frame), x, ...
                'AbsTol', 1e-12);
        end

        function earlyCorridorSwitchKeepsExplicitIncomingFrame(testCase)
            [model, corridors] = testCase.fixture();
            z = [model.CumulativeLength(2); 0.2; -1; 0; 0; 0];
            [A, b, corridor] = corridors.getPathConstraintsAtS(z(1), 1);
            testCase.verifyEqual(corridor, 2);
            x = model.pathToCartesian(z, 1);
            [Axy, bxy] = corridors.getCartesianConstraints(corridor);
            testCase.verifyEqual(A*z-b, Axy*x(1:2)-bxy, 'AbsTol', 1e-11);
        end

        function finalFrameExtrapolationMatchesCartesianFaces(testCase)
            [model, corridors] = testCase.fixture();
            z = [model.TotalLength+2; -0.2; 0.7; 1; 0; 0];
            [A, b, corridor] = corridors.getPathConstraintsAtS(z(1), 3);
            testCase.verifyEqual(corridor, 3);
            x = model.pathToCartesian(z, 3);
            [Axy, bxy] = corridors.getCartesianConstraints(corridor);
            testCase.verifyEqual(A*z-b, Axy*x(1:2)-bxy, 'AbsTol', 1e-11);
        end

        function controllerUsesMeasuredFrameForCorridorCheck(testCase)
            [model, corridors, config] = testCase.fixture();
            controller = NominalMPC(model, corridors, config);
            x = [15; 3; 0; 1; 0; 0];
            [~, info] = controller.computeControl(x);
            testCase.verifyEqual(info.FrameSegmentIndex, 2);
            testCase.verifyLessThan(info.PathState(3), -4.49);
            [Axy, bxy] = corridors.getCartesianConstraints(info.SegmentIndex);
            testCase.verifyEqual(info.InsideCorridor, ...
                all(Axy*x(1:2) <= bxy+1e-8));
            testCase.verifyEqual( ...
                model.pathToCartesian(info.PathState, info.FrameSegmentIndex), ...
                x, 'AbsTol', 1e-12);
        end

        function goalReturnIncludesFrame(testCase)
            [model, corridors, config] = testCase.fixture();
            controller = NominalMPC(model, corridors, config);
            model.updateFromCartesian([10; 13; 0; 0; 0; 0]);
            x = [model.Path(end, :)'; 0; 0; 0; 0];
            [u, info] = controller.computeControl(x);
            testCase.verifyTrue(info.GoalReached);
            testCase.verifyEqual(info.FrameSegmentIndex, 3);
            testCase.verifyEqual(u, zeros(2, 1));
            testCase.verifyEqual( ...
                model.pathToCartesian(info.PathState, info.FrameSegmentIndex), ...
                x, 'AbsTol', 1e-12);
        end

        function currentSimulationDetectsErrorAtFirstTurn(testCase)
            [model, corridors, config, parameters] = testCase.fixture();
            controller = NominalMPC(model, corridors, config);
            plant = AUVVehicle(parameters);
            foundTurn = false;
            for k = 1:120
                [u, info] = controller.computeControl(plant.State);
                reconstructed = model.pathToCartesian( ...
                    info.PathState, info.FrameSegmentIndex);
                testCase.verifyEqual(reconstructed(1:2), plant.State(1:2), ...
                    'AbsTol', 1e-11);
                if info.FrameSegmentIndex == 2 && plant.State(1) > 10.5
                    testCase.verifyLessThan(info.PathState(3), -0.01);
                    testCase.verifyGreaterThan(u(2), 0);
                    foundTurn = true;
                    break
                end
                plant.step(u);
            end
            testCase.verifyTrue(foundTurn, ...
                'The current simulation must detect lateral error at its first turn.');
        end
    end

    methods (Static, Access = private)
        function [model, corridors, config, p] = fixture()
            p = struct('mx', 30, 'my', 40, 'Iz', 8, ...
                'dx', 12, 'dy', 20, 'dr', 5, 'SampleTime', 0.1, ...
                'InitialState', [2; 3; 0; 0; 0; 0], ...
                'StateCovariance', zeros(6));
            path = [2 3; 10.5 3; 10 13; 32 6];
            model = AUVPathModel(path, p);
            map = binaryOccupancyMap(40, 15, 1);
            setOccupancy(map, [8 9; 14 9; 18 7; 26 11; 30 3], ones(5, 1));
            corridors = LiuSafeCorridor(0.4, 3);
            corridors.generate(path, map);
            config = struct('N', 10, ...
                'Q', diag([0.5; 10; 500; 1; 1; 2]), ...
                'R', diag([0.001; 0.005]), 'ReferenceSpeed', 3, ...
                'Umin', [-200; -20], 'Umax', [200; 20], ...
                'GoalTolerance', 0.03);
        end
    end
end
