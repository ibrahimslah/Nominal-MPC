# Nominal MPC for AUV Path Following

MATLAB implementation of nominal model predictive control (MPC) for a planar autonomous underwater vehicle (AUV). The controller follows a supplied piecewise-linear path through safe corridors generated from a binary occupancy map, while respecting thrust and yaw-moment bounds.

The simulation applies MPC inputs to a nonlinear vehicle model. The controller linearizes its path-coordinate prediction model at each control step and solves a quadratic program using `quadprog`.

## Requirements

- MATLAB.
- Navigation Toolbox for `binaryOccupancyMap`.
- Optimization Toolbox for `quadprog` and `optimoptions`.

The standalone safe-corridor example requires Navigation Toolbox but does not require Optimization Toolbox. Tests use MATLAB's `matlab.unittest` framework. A minimum MATLAB release has not been established.

## Quick start

Set MATLAB's current folder to this project folder, then run:

```matlab
runNominalMPCSimulation
```

The script adds the project and `SafeCorridor` folders to the MATLAB path, creates an occupancy map and waypoint path, generates corridors, and runs the closed-loop simulation. It begins with `clear`, `clc`, and `close all`.

Three figures show:

1. The occupancy map, path, safe corridors, and AUV trajectory.
2. Path progress, lateral error, and surge velocity over time.
3. Applied thrust and yaw moment over time.

To explore corridor generation independently:

```matlab
run(fullfile('SafeCorridor', 'LiuSafeCorridorExample.m'))
```

This example demonstrates Cartesian and path-coordinate inequality queries, online corridor selection, and visualization on a map with a shifted world origin.

## Project structure

```text
Nominal MPC/
├── runNominalMPCSimulation.m        Closed-loop example and plots
├── AUVVehicle.m                     Nonlinear Cartesian plant
├── AUVPathModel.m                   Path coordinates and prediction model
├── @NominalMPC/
│   ├── NominalMPC.m                 Stateful controller and QP solve
│   ├── liftedAffine.m               Horizon dynamics assembly
│   ├── buildNominalMPCCost.m        Quadratic tracking cost
│   └── buildPathMPCConstraints.m    Horizon corridor inequalities
├── SafeCorridor/
│   ├── LiuSafeCorridor.m            Corridor generation and queries
│   ├── LiuSafeCorridorExample.m     Standalone corridor example
│   └── tests/                      Corridor geometry and API tests
└── tests/                          Path-model and frame-integration tests
```

`AUVVehicle` integrates the nonlinear plant with `ode45` and can add Gaussian process noise after each step. `AUVPathModel` manages path geometry, coordinate transformations, and nonlinear and linearized discrete predictions. `LiuSafeCorridor` generates and caches one corridor per path segment. `NominalMPC` manages control memory and returns an input; the caller advances the plant and records results.

## States and inputs

The plant state is `xCartesian = [px; py; yaw; vx; vy; r]`. Positions are in world coordinates; velocities are in the vehicle body frame.

The controller uses `xPath = [s; ePsi; eY; vx; vy; r]`:

| State | Meaning | Unit |
| --- | --- | --- |
| `s` | Global cumulative progress along the path | m |
| `ePsi` | Heading error relative to the segment frame | rad |
| `eY` | Signed lateral displacement from the segment centerline | m |
| `vx` | Surge velocity | m/s |
| `vy` | Sway velocity | m/s |
| `r` | Yaw rate | rad/s |

The control input is `u = [Fx; tau]`, with surge thrust `Fx` in N and yaw moment `tau` in N·m.

## Controller setup and use

The following example is self-contained when run from the project folder:

```matlab
addpath(pwd);
addpath(fullfile(pwd, 'SafeCorridor'));

p = struct('mx', 30, 'my', 40, 'Iz', 8, ...
    'dx', 12, 'dy', 20, 'dr', 5, 'SampleTime', 0.1, ...
    'InitialState', [2; 3; 0; 0; 0; 0], ...
    'StateCovariance', zeros(6));
path = [2 3; 10 3; 20 6; 32 6];
map = binaryOccupancyMap(40, 15, 1);

plant = AUVVehicle(p);
pathModel = AUVPathModel(path, p);
corridorGenerator = LiuSafeCorridor(0.4, 3.0);
corridorGenerator.generate(path, map);

config = struct('N', 100, ...
    'Q', diag([0.5; 10; 50; 10; 1; 2]), ...
    'R', diag([0.01; 0.05]), 'ReferenceSpeed', 3, ...
    'Umin', [-200; -20], 'Umax', [200; 20], ...
    'GoalTolerance', 0.03);
controller = NominalMPC(pathModel, corridorGenerator, config);

[uk, info] = controller.computeControl(plant.State);
if ~info.GoalReached
    xNext = plant.step(uk);
end
```

This setup uses an empty map; the main simulation adds obstacles. Generate corridors before constructing the controller, using exactly the same path as the path model.

`computeControl` returns a two-element input and an `info` structure containing the measured `PathState`, frame and corridor indices, corridor-membership check, reference trajectory, solver status, cost, and fallback/goal flags. Successful solves also provide the stacked `OptimalSequence` and `PredictedState`. Predictions contain the current state and `N` future states.

`info.FrameSegmentIndex` identifies the coordinate frame used by the measured path state and every predicted state in that solve. `info.SegmentIndex` identifies the selected corridor. These indices can differ during corridor transitions. When converting a logged path state back to Cartesian coordinates, preserve its frame:

```matlab
xCartesian = pathModel.pathToCartesian( ...
    info.PathState, info.FrameSegmentIndex);
```

`controller.reset()` clears the previous input and stored control sequence, leaving the path model and generated corridors intact.

## Configuration

Edit `runNominalMPCSimulation.m` to change the scenario and controller settings.

| Setting | Main simulation value | Purpose |
| --- | --- | --- |
| `path` | Four waypoints from `[2, 3]` to `[32, 6]` | Supplied path in world coordinates |
| Map dimensions and resolution | 40 m × 15 m, 1 cell/m | Workspace and obstacle discretization |
| `obstacles` | Five occupied world-coordinate points | Occupied cells in the map |
| `p.SampleTime` | 0.1 s | Plant and controller sample interval |
| `p.InitialState` | `[2; 3; 0; 0; 0; 0]` | Initial Cartesian plant state |
| `p.StateCovariance` | Zero 6 × 6 matrix | Per-step additive plant-noise covariance |
| `robotRadius` | 0.4 m | Vehicle clearance parameter for corridors |
| `boundingBoxMargin` | 3.0 m | Local corridor bounding-box margin |
| `config.N` | 100 | Prediction horizon in control steps |
| `config.Q` | `diag([0.5; 10; 50; 10; 1; 2])` | Weights in path-state order |
| `config.R` | `diag([0.01; 0.05])` | Weights in input order |
| `config.ReferenceSpeed` | 3 m/s | Reference surge speed and progress rate |
| `config.Umin`, `config.Umax` | `[-200; -20]`, `[200; 20]` | Thrust and moment bounds |
| `config.GoalTolerance` | 0.03 m | Progress threshold for goal detection |
| `maxSimulationSteps` | 60000 | Maximum closed-loop iterations |

Vehicle parameters `mx`, `my`, and `Iz` set inertial terms; `dx`, `dy`, and `dr` set damping terms.

Required controller fields are `N`, `Q`, `R`, `ReferenceSpeed`, `Umin`, and `Umax`. `Q` must be symmetric positive semidefinite and `R` symmetric positive definite. Input bounds must be ordered and contain zero. Optional `GoalTolerance` defaults to `0.03`. Optional `SolverOptions` must be `quadprog` options using the `active-set` algorithm; the default also disables solver display.

## Simulation outputs and behavior

The main script retains these arrays in the MATLAB workspace:

| Variable | Contents |
| --- | --- |
| `xCartesianHistory` | 6 × state-sample count Cartesian states |
| `xPathHistory` | 6 × state-sample count path states |
| `uHistory` | 2 × applied-control count inputs |
| `segmentHistory` | Corridor index for each state sample |
| `frameSegmentHistory` | Coordinate-frame index for each state sample |
| `stateTime`, `controlTime` | Time vectors for states and inputs |

The script also retains generated `corridors` and the plant, path model, and controller objects. It creates figures but does not automatically save plots or histories to files.

At each solve, the controller applies the selected corridor inequalities to all future prediction steps in the measured segment frame, with input bounds over the horizon. The current measured state is checked for corridor membership separately. The simulation warns if it lies outside the selected corridor.

When `quadprog` returns a nonpositive exit flag, the controller uses the next input from its stored, shifted solution. Without a stored solution, it returns zero input. The simulation warns and continues. This fallback is not a new feasibility check and does not guarantee that the actual vehicle remains inside a corridor.

Goal detection uses `s >= TotalLength - GoalTolerance`. The controller returns zero input and clears its control memory; the simulation stops before applying that input. This is a progress criterion, rather than a terminal position-and-velocity check.

The example sets process noise to zero. Although the plant supports nonzero Gaussian noise, the controller is nominal: it does not implement covariance steering, chance constraints, or probabilistic collision guarantees. Paths are supplied by the caller; this example does not search for a route through the map.

## Running tests

From the project folder in MATLAB:

```matlab
addpath(pwd);
addpath(fullfile(pwd, 'SafeCorridor'));
results = [runtests('tests'); ...
           runtests(fullfile('SafeCorridor', 'tests'))];
disp(results);
assertSuccess(results);
```

The tests cover path-coordinate transformations and segment-frame handling, controller/frame integration, corridor geometry, Cartesian/path inequality consistency, online corridor transitions, and the corridor public API. Run them in your MATLAB environment to verify compatibility and behavior.
