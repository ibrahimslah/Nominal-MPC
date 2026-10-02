%% Nominal MPC simulation for the nonlinear AUV plant.

clear;
clc;
close all;

% Allow the example to run from a different working directory.
exampleDirectory = fileparts(mfilename('fullpath'));
addpath(exampleDirectory);
addpath(fullfile(exampleDirectory, 'SafeCorridor'));

%% Vehicle parameters

p.mx = 30;
p.my = 40;
p.Iz = 8;
p.dx = 12;
p.dy = 20;
p.dr = 5;
p.SampleTime = 0.1;
p.InitialState = [2; 3; 0; 0; 0; 0];
p.StateCovariance = diag([0; 0; 0; 0; 0; 0]);

Ts = p.SampleTime;

%% Path and occupancy map

path = [
     2    3
    10    3
    20    6
    32    6
];

mapWidth = 40;
mapHeight = 15;
resolution = 1;
map = binaryOccupancyMap(mapWidth, mapHeight, resolution);

obstacles = [
     8    9
    14    9
    18    7
    26   11
    30    3
];
setOccupancy(map, obstacles, ones(size(obstacles, 1), 1));

%% Plant, path model, and safe corridors

plant = AUVVehicle(p);
pathModel = AUVPathModel(path, p);

robotRadius = 0.4;
boundingBoxMargin = 3.0;
corridorGenerator = LiuSafeCorridor(robotRadius, boundingBoxMargin);
corridors = corridorGenerator.generate(path, map);

%% MPC configuration

nx = 6;
nu = 2;
N = 100;
vRef =3;

% Path state: [s; ePsi; eY; vx; vy; r]. Input: [Fx; tau].
Q = diag([.5; 10; 50; 10; 1; 2]);
R = diag([0.01; .05]);

FxMin = -200;
FxMax = 200;
tauMin = -20;
tauMax = 20;
umin = [FxMin; tauMin];
umax = [FxMax; tauMax];

config = struct( ...
    'N', N, ...
    'Q', Q, ...
    'R', R, ...
    'ReferenceSpeed', vRef, ...
    'Umin', umin, ...
    'Umax', umax, ...
    'GoalTolerance', 0.03, ...
    'SolverOptions', optimoptions('quadprog', ...
        'Algorithm', 'active-set', 'Display', 'off'));

controller = NominalMPC(pathModel, corridorGenerator, config);

%% Simulation and logging

maxSimulationSteps = 60000;
xCartesianHistory = zeros(nx, maxSimulationSteps + 1);
xPathHistory = zeros(nx, maxSimulationSteps + 1);
uHistory = zeros(nu, maxSimulationSteps);
segmentHistory = zeros(1, maxSimulationSteps + 1);
% Keep the coordinate frame separately from the selected corridor index.
frameSegmentHistory = zeros(1, maxSimulationSteps + 1);
xCartesianHistory(:, 1) = plant.State;

numberOfSteps = 0;
numberOfStateSamples = 0;

for k = 1:maxSimulationSteps
    [uk, info] = controller.computeControl(plant.State);

    % This measurement corresponds to the state before input k is applied.
    xPathHistory(:, k) = info.PathState;
    segmentHistory(k) = info.SegmentIndex;
    frameSegmentHistory(k) = info.FrameSegmentIndex;
    numberOfStateSamples = k;

    if info.GoalReached
        fprintf('Goal reached after %d control steps.\n', numberOfSteps);
        break;
    end

    if ~info.InsideCorridor
        warning('Current state outside corridor at MPC step %d.', k);
    end

    if info.UsedFallback
        warning('QP failed at step %d. exitflag = %d', k, info.ExitFlag);
    end

    xNext = plant.step(uk);
    uHistory(:, k) = uk;
    xCartesianHistory(:, k + 1) = xNext;
    numberOfSteps = k;
end

% If the step limit was reached, also log its final post-step state.
if numberOfStateSamples == numberOfSteps
    [finalPathState, finalFrameSegment] = ...
        pathModel.updateFromCartesian(plant.State);
    xPathHistory(:, numberOfSteps + 1) = finalPathState;
    frameSegmentHistory(numberOfSteps + 1) = finalFrameSegment;
    [~, ~, segmentHistory(numberOfSteps + 1)] = ...
        corridorGenerator.getPathConstraintsAtS(finalPathState(1), finalFrameSegment);
    numberOfStateSamples = numberOfSteps + 1;
end
%%
xCartesianHistory = xCartesianHistory(:, 1:numberOfSteps + 1);
xPathHistory = xPathHistory(:, 1:numberOfStateSamples);
uHistory = uHistory(:, 1:numberOfSteps);
segmentHistory = segmentHistory(1:numberOfStateSamples);
frameSegmentHistory = frameSegmentHistory(1:numberOfStateSamples);

%% Map, corridor, path, and vehicle trajectory

ax = corridorGenerator.plot();
hold(ax, 'on');
plot(ax, xCartesianHistory(1, :), xCartesianHistory(2, :), ...
    'b-', 'LineWidth', 2, 'DisplayName', 'AUV trajectory');
plot(ax, xCartesianHistory(1, 1), xCartesianHistory(2, 1), ...
    'go', 'MarkerSize', 8, 'LineWidth', 2, ...
    'DisplayName', 'Start position');
plot(ax, xCartesianHistory(1, end), xCartesianHistory(2, end), ...
    'bo', 'MarkerSize', 8, 'LineWidth', 2, ...
    'DisplayName', 'Final position');
legend(ax, 'Location', 'best');
hold(ax, 'off');

%% Path-coordinate states

stateTime = (0:numberOfStateSamples - 1) * Ts;
figure;

subplot(3, 1, 1);
plot(stateTime, xPathHistory(1, :), 'LineWidth', 1.5);
ylabel('s [m]');
grid on;

subplot(3, 1, 2);
plot(stateTime, xPathHistory(3, :), 'LineWidth', 1.5);
ylabel('e_Y [m]');
grid on;

subplot(3, 1, 3);
plot(stateTime, xPathHistory(4, :), 'LineWidth', 1.5);
yline(vRef, '--');
ylabel('v_x [m/s]');
xlabel('Time [s]');
grid on;

%% Inputs

controlTime = (0:numberOfSteps - 1) * Ts;
figure;

subplot(2, 1, 1);
stairs(controlTime, uHistory(1, :), 'LineWidth', 1.5);
ylabel('F_x [N]');
grid on;

subplot(2, 1, 2);
stairs(controlTime, uHistory(2, :), 'LineWidth', 1.5);
ylabel('\tau [Nm]');
xlabel('Time [s]');
grid on;
