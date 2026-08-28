function C = solve_nutrient(rho, cfg)
%SOLVE_NUTRIENT Steady-state reaction-diffusion field on a coarse grid.
%
%   C = SOLVE_NUTRIENT(rho, cfg) solves
%
%       D * laplacian(C) = k * rho * C      <=>      L^2 laplacian(C) = rho*C
%
%   with L = sqrt(D/k) the diffusion length (fine sites), rho in [0,1] the
%   occupied-cell fraction per coarse cell, and Dirichlet boundary C = C0 on
%   the frame edge.
%
%   The equation is solved directly as a sparse linear system. Explicit
%   relaxation to a diffusive steady state needs ~N^2 iterations and is not
%   feasible at 1024^2; the sparse solve is O((nc^2)^1.5) and lands in tens
%   of milliseconds.
%
%   Discretisation (5-point Laplacian, spacing hc = N/nc, lambda = L^2/hc^2):
%
%       lambda * (C_up + C_down + C_left + C_right - 4*C) = rho*C

n  = size(rho, 1);
hc = cfg.sim.N / cfg.sim.nc;
lam = (cfg.sim.L / hc)^2;
C0 = cfg.sim.C0;

% Interior mask: everything but the one-cell Dirichlet ring.
mask = true(n);
mask(1, :) = false; mask(end, :) = false;
mask(:, 1) = false; mask(:, end) = false;

N2 = n * n;
intLin = find(mask);
bound  = find(~mask);

[ii, jj] = ind2sub([n n], intLin);
rhoFlat  = rho(intLin);

up    = sub2ind([n n], ii - 1, jj);
down  = sub2ind([n n], ii + 1, jj);
left  = sub2ind([n n], ii,     jj - 1);
right = sub2ind([n n], ii,     jj + 1);
nbrs  = [up down left right];            % M x 4

nbrIsInt = mask(nbrs);                    % M x 4 logical

% Diagonal (centre) coefficients, one per interior row.
diagVal = -(4 * lam + rhoFlat);

% Off-diagonal lam for each interior neighbour. Boundary neighbours are
% Dirichlet (C = C0) and move to the RHS.
rows = repmat(intLin, 1, 4);
cols = nbrs;
vals = lam * ones(size(nbrs));
sel  = nbrIsInt(:);

A = sparse([intLin; rows(sel)], [intLin; cols(sel)], ...
           [diagVal; vals(sel)], N2, N2);

% Dirichlet rows: C = C0 on the ring.
A = A + sparse(bound, bound, 1, N2, N2);

b = zeros(N2, 1);
b(bound)  = C0;
b(intLin) = -lam * C0 * sum(~nbrIsInt, 2);   % boundary-neighbour contribution

Csol = A \ b;
C = reshape(Csol, n, n);
end
