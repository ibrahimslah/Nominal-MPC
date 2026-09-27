function [Alift, Blift, rlift] = liftedAffine(A, B, r, N)
%LIFTEDAFFINE Build lifted affine prediction model.
%
%   [Alift, Blift, rlift] = liftedAffine(A, B, r, N)
%
% Builds the lifted prediction model
%
%       X = Alift*x0 + Blift*U + rlift
%
% for the affine discrete-time system
%
%       x(k+1) = A(k)*x(k) + B(k)*u(k) + r(k)
%
% using the stacking convention
%
%       X = [x0;
%            x1;
%            ...
%            xN]
%
%       U = [u0;
%            u1;
%            ...
%            u(N-1)]
%
% This is the same state-stacking convention used in the
% covariance-steering formulation:
%
%       X = A*x0 + B*U + D*W + R
%
% The disturbance term D*W is not included here.
%
%
% INPUTS
% -------------------------------------------------------------------------
% A : State-transition matrix/matrices.
%
%     LTI:
%         A is nx-by-nx
%
%     LTV:
%         A(:,:,k) is nx-by-nx, k = 1,...,N
%
% B : Input matrix/matrices.
%
%     LTI:
%         B is nx-by-nu
%
%     LTV:
%         B(:,:,k) is nx-by-nu, k = 1,...,N
%
% r : Affine term.
%
%     LTI:
%         r is nx-by-1
%
%     LTV:
%         r(:,k) is nx-by-1, k = 1,...,N
%
% N : Prediction horizon
%
%
% OUTPUTS
% -------------------------------------------------------------------------
% Alift : ((N+1)*nx)-by-nx
%
% Blift : ((N+1)*nx)-by-(N*nu)
%
% rlift : ((N+1)*nx)-by-1
%
%
% The resulting prediction is
%
%       X = Alift*x0 + Blift*U + rlift
%
% where X contains x0,...,xN.


%% Dimensions

nx = size(A,1);
nu = size(B,2);


%% Allocate lifted matrices

Alift = zeros((N+1)*nx, nx);

Blift = zeros((N+1)*nx, N*nu);

rlift = zeros((N+1)*nx, 1);


%% Initial state
%
% x0 = I*x0

Alift(1:nx,:) = eye(nx);

% Blift first block row = 0
% rlift first block     = 0


%% Current propagated mappings
%
% At prediction step i:
%
%   xi = Phi*x0 + Gamma*U + eta

Phi   = eye(nx);
Gamma = zeros(nx, N*nu);
eta   = zeros(nx,1);


%% Build prediction model recursively

for i = 1:N

    % -------------------------------------------------------------
    % Select model at prediction step i
    % -------------------------------------------------------------

    if ismatrix(A)
        Ai = A;
    else
        Ai = A(:,:,i);
    end

    if ismatrix(B)
        Bi = B;
    else
        Bi = B(:,:,i);
    end

    if size(r,2) == 1
        ri = r;
    else
        ri = r(:,i);
    end


    % -------------------------------------------------------------
    % Propagate state mapping
    %
    % x_i = Ai*x_{i-1} + Bi*u_{i-1} + ri
    % -------------------------------------------------------------

    Phi = Ai * Phi;

    Gamma = Ai * Gamma;

    inputColumns = ...
        (i-1)*nu + (1:nu);

    Gamma(:,inputColumns) = ...
        Gamma(:,inputColumns) + Bi;

    eta = Ai * eta + ri;


    % -------------------------------------------------------------
    % Store lifted blocks
    % -------------------------------------------------------------

    stateRows = ...
        i*nx + (1:nx);

    Alift(stateRows,:) = Phi;

    Blift(stateRows,:) = Gamma;

    rlift(stateRows,:) = eta;

end

end