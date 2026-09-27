function [Aineq, bineq] = buildPathMPCConstraints( ...
    xk, Alift, Blift, rlift, A_path, b_path, N)
%BUILDPATHMPCCONSTRAINTS
%
% Lifted dynamics:
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
%
% Path constraints:
%
%   x0 : no constraint
%
%        0*x0 <= 0
%
%   xi : A_path*xi <= b_path
%        i = 1,...,N
%
%
% Result:
%
%   Aineq*U <= bineq


    %% Dimensions

    nx = size(Alift,2);

    nc = size(A_path,1);


    %% ================================================================
    % Build lifted state-constraint matrix
    %
    % Ac =
    %
    % [ 0          0          ...       0
    %   0        A_path       ...       0
    %   ...
    %   0          0          ...    A_path ]
    %
    % corresponding to
    %
    % X = [x0; x1; ...; xN]
    % =================================================================

    A0 = zeros(nc,nx);

    Abar = blkdiag( ...
        A0, ...
        kron(speye(N),A_path) );


    %% ================================================================
    % Lifted right-hand side
    %
    % bbar =
    %
    % [0;
    %  b_path;
    %  ...
    %  b_path]
    % =================================================================

    bbar = [
        zeros(nc,1)
        repmat(b_path,N,1)
    ];


    %% ================================================================
    % Substitute lifted dynamics
    %
    % Abar * (Alift*xk + Blift*U + rlift) <= bbar
    %
    % Therefore:
    %
    % Abar*Blift*U <=
    %
    % bbar - Abar*(Alift*xk + rlift)
    % =================================================================

    Aineq = ...
        Abar * Blift;

    bineq = ...
        bbar ...
        - Abar*( ...
            Alift*xk ...
            + rlift );

end