%% LiuSafeCorridor public API example
% This example requires Navigation Toolbox for binaryOccupancyMap.
% Optimization Toolbox is not required.

exampleFolder = fileparts(mfilename('fullpath'));
addpath(exampleFolder);

%% Shifted occupancy map, obstacles, and piecewise-linear path
map = binaryOccupancyMap(false(120, 160), 10);
map.GridLocationInWorld = [-20 30];

obstaclePoints = [ ...
    -16 36
    -12 34
     -8 39
     -5 35];
setOccupancy(map, obstaclePoints, ...
    true(size(obstaclePoints, 1), 1));

path = [ ...
    -18 32
    -13 32
     -9 36
     -5 38];

%% Generate and cache one safe corridor per path segment
generator = LiuSafeCorridor(0.25, 2.5);
corridors = generator.generate(path, map);

%% Query one segment directly
segmentIndex = 1;
[Axy, bxy] = generator.getCartesianConstraints(segmentIndex);
[Apath, bpath] = generator.getPathConstraints(segmentIndex);

% Axy*[x;y] <= bxy.
% Apath*[s;ePsi;eY;vx;vy;r] <= bpath.
fprintf('Generated %d safe corridors.\n', numel(corridors));
fprintf(['Segment %d has %d Cartesian inequalities (%d columns) and ' ...
    '%d path inequalities (%d columns).\n'], ...
    segmentIndex, numel(bxy), size(Axy, 2), ...
    numel(bpath), size(Apath, 2));

%% Select constraints online using global cumulative progress
sQuery = mean(corridors(2).SLimits);
[ApathAtS, bpathAtS, selectedSegment] = ...
    generator.getPathConstraintsAtS(sQuery);
fprintf(['s = %.3f m selects segment %d with %d inequalities ' ...
    'and %d state columns.\n'], ...
    sQuery, selectedSegment, numel(bpathAtS), size(ApathAtS, 2));

%% Visualize the cached map, path, and corridors
ax = generator.plot();
title(ax, 'Liu safe-corridor public API example');
drawnow;
