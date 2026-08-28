function S = step_migrate(S)
%STEP_MIGRATE Random walk gated by adhesion and steered by chemotaxis.
%
%   Migration probability is
%
%       p = P_move * (1 - A_i * (n_occ / 4))
%
%   High adhesion A_i keeps well-attached cells in the aggregate; low A_i
%   lets cells detach and wander, forming satellite clusters. Direction is a
%   uniform random empty neighbour, optionally biased up the nutrient
%   gradient by chi. Note: the most distinctive migration phenotype --
%   satellite detachment -- is exactly what Part 1's QC pipeline strips or
%   rejects, so +bridge records the stripped-area fraction separately.

cfg = S.cfg;
N   = S.N;
occ = S.occ;
age = S.age;
clone = S.clone;

P_move = cfg.sim.P_move;
A_i    = cfg.sim.A_i;
conn   = cfg.sim.connectivity;
diagW  = cfg.sim.diagWeight;
chi    = cfg.sim.chi;
Cfine  = S.Cfine;

kern = [0 1 0; 1 0 1; 0 1 0];
ncount = conv2(double(occ), kern, 'same');

empty = ~occ;
if conn == 4
    se = strel('disk', 1, 0);
else
    se = strel('square', 3);
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

    c = ncount(lin);
    p = P_move * (1 - A_i * (c / 4));
    if p <= 0, continue; end

    if rand() < p
        [nbr, ok] = ca.sample_empty_neighbor(lin, occ, N, conn, diagW, Cfine, chi);
        if ok
            occ(lin)   = false;
            occ(nbr)   = true;
            age(nbr)   = age(lin);
            age(lin)   = 0;
            clone(nbr) = clone(lin);
            clone(lin) = 0;
        end
    end
end

S.occ = occ;
S.age = age;
S.clone = clone;
end
