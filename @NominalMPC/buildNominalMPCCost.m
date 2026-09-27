function [H, F] = buildNominalMPCCost(~, xk, Alift, Blift, rlift, Xref, Q, R, N)
%BUILDNOMINALMPCCOST
%
% Build the quadratic cost matrices for nominal MPC.
%
% Lifted prediction:
%
%   X = Alift*xk + Blift*U + rlift
%
% where
%
%   X = [x0;
%        x1;
%        ...
%        xN]
%
%   U = [u0;
%        u1;
%        ...
%        u_{N-1}]
%
%
% MPC cost:
%
%   J = (X - Xref)'*Qbar*(X - Xref)
%       + U'*Rbar*U
%
%
% quadprog form:
%
%   min  0.5*U'*H*U + F'*U
%
%
% INPUTS
% -------------------------------------------------------------------------
% xk      : nx x 1
%           Current path-coordinate state
%
% Alift   : (N+1)*nx x nx
%
% Blift   : (N+1)*nx x N*nu
%
% rlift   : (N+1)*nx x 1
%
% Xref    : (N+1)*nx x 1
%
% Q       : nx x nx
%           State tracking weight
%
% R       : nu x nu
%           Input weight
%
% N       : prediction horizon
%
%
% OUTPUTS
% -------------------------------------------------------------------------
% H       : N*nu x N*nu
%
% F       : N*nu x 1
%
% such that
%
%   J_QP = 0.5*U'*H*U + F'*U
%
% differs from the original cost only by a constant term.


    %% ================================================================
    % Dimensions
    % =================================================================

    nx = size(Q,1);
    nu = size(R,1);


    %% ================================================================
    % Lifted state weighting matrix
    %
    % X = [x0; x1; ...; xN]
    %
    % x0 is the measured/current state and cannot be changed by U.
    %
    % Therefore:
    %
    % Qbar =
    %
    % [ 0   0   ...   0
    %   0   Q   ...   0
    %   ...
    %   0   0   ...   Q ]
    %
    % This preserves the same lifted dimension while giving x0
    % zero optimization weight.
    % =================================================================

    Qbar = blkdiag( ...
        zeros(nx,nx), ...
        kron(speye(N),Q) );


    %% ================================================================
    % Lifted input weighting matrix
    %
    % U = [u0; ...; u_{N-1}]
    % =================================================================

    Rbar = ...
        kron(speye(N),R);


    %% ================================================================
    % Free state trajectory
    %
    % This is the predicted trajectory when U = 0:
    %
    %   Xfree = Alift*xk + rlift
    % =================================================================

    Xfree = ...
        Alift*xk ...
        + rlift;


    %% ================================================================
    % Tracking error independent of optimization variable U
    %
    % X - Xref =
    %
    %   Blift*U + error
    % =================================================================

    error = ...
        Xfree ...
        - Xref;


    %% ================================================================
    % Expand cost
    %
    % J =
    %
    % (Blift*U + error)' Qbar (Blift*U + error)
    % + U' Rbar U
    %
    %
    % =
    %
    % U' (Blift'*Qbar*Blift + Rbar) U
    %
    % + 2 error'*Qbar*Blift U
    %
    % + constant
    %
    %
    % quadprog uses:
    %
    %   0.5 U'*H*U + F'*U
    %
    % Therefore:
    %
    % H = 2*(Blift'*Qbar*Blift + Rbar)
    %
    % F = 2*Blift'*Qbar*error
    % =================================================================

    H = 2*( ...
        Blift.'*Qbar*Blift ...
        + Rbar );


    F = 2*( ...
        Blift.'*Qbar*error );


    %% ================================================================
    % Numerical symmetry
    %
    % Theoretically H is symmetric.
    % This removes small floating-point asymmetry.
    % =================================================================

    H = (H + H.')/2;


    %% ================================================================
    % Dimension checks
    % =================================================================

    assert( ...
        size(Alift,1) == (N+1)*nx, ...
        'Alift has inconsistent dimensions.');

    assert( ...
        size(Blift,1) == (N+1)*nx, ...
        'Blift has inconsistent dimensions.');

    assert( ...
        length(rlift) == (N+1)*nx, ...
        'rlift has inconsistent dimensions.');

    assert( ...
        length(Xref) == (N+1)*nx, ...
        'Xref must contain [xref0; ...; xrefN].');

    assert( ...
        size(H,1) == N*nu, ...
        'H has inconsistent dimensions.');

end
