function S = step_proliferate(S)
%STEP_PROLIFERATE Division gated by maturity, space, contact inhibition, nutrient.
%
%   A cell divides only if (a) age >= T_cell, (b) an empty neighbour exists,
%   (c) contact inhibition beta does not block, and (d) local nutrient is
%   adequate. Division probability is
%
%       p = alpha_clone * (1 - beta * (n_occ / 4)) * min(1, C / C_q)
%
%   where n_occ in 0..4 counts occupied ORTHOGONAL neighbours. Benign beta~1
%   blocks division in crowded sites; malignant beta~0 lets a surrounded cell
%   push into its one remaining gap, which is the mechanistic source of the
%   infiltrative / lobular margin.

cfg = S.cfg;
N   = S.N;
occ = S.occ;
age = S.age;
clone = S.clone;

T_cell = cfg.sim.T_cell;
beta   = cfg.sim.beta;
conn   = cfg.sim.connectivity;
diagW  = cfg.sim.diagWeight;
Cq     = cfg.sim.C_q;
Cfine  = S.Cfine;
alphaByClone = S.alphaByClone;

% Occupied orthogonal-neighbour count (0..4) at every occupied site.
kern = [0 1 0; 1 0 1; 0 1 0];
ncount = conv2(double(occ), kern, 'same');

% Frontier = occupied cells with at least one empty neighbour.
empty = ~occ;
if conn == 4
    se = strel('disk', 1, 0);            % plus-shaped (4-connectivity)
else
    se = strel('square', 3);             % full 3x3 (8-connectivity)
end
frontierMask = occ & imdilate(empty, se);
frontier = find(frontierMask);

if isempty(frontier)
    S.occ = occ; S.age = age; S.clone = clone;
    return
end

order = randperm(numel(frontier));

for k = 1:numel(order)
    lin = frontier(order(k));
    if ~occ(lin), continue; end
    if age(lin) < T_cell, continue; end

    c = ncount(lin);                       % occupied orthogonal neighbours

    a = alphaByClone(clone(lin));
    p = a * (1 - beta * (c / 4));
    if p <= 0, continue; end

    % Nutrient gating (soft quiescence below C_q).
    if cfg.sim.nutrientCoupling && isfinite(Cq) && Cq > 0
        p = p * min(1, Cfine(lin) / Cq);
    end

    if rand() < p
        [nbr, ok] = ca.sample_empty_neighbor(lin, occ, N, conn, diagW, Cfine, 0);
        if ok
            occ(nbr)   = true;
            age(nbr)   = 0;
            age(lin)   = 0;               % T_cell is a MINIMUM interval
            clone(nbr) = clone(lin);
            S.totalDivisions = S.totalDivisions + 1;
        end
    end
end

S.occ = occ;
S.age = age;
S.clone = clone;
end
