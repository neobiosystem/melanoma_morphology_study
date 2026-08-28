function S = step_proliferate(S)
%STEP_PROLIFERATE Division + the two new Part 3 mechanisms.
%
%   Part 2's rule (division gated by maturity, space, contact inhibition,
%   nutrient) with two extensions that fill the "surrounded cell cannot
%   divide" gap handed over by Part 2:
%
%     (a) push displacement -- when cfg.sim.displace is true, cells one layer
%         in from the frontier (surrounded, but adjacent to a frontier cell)
%         are also eligible: on passing the division test they push a
%         neighbour outward via MECH.PUSH_DISPLACE instead of dropping the
%         division. This lets a second shell divide (bulk growth), not just
%         the Eden surface.
%     (b) growth-rate field -- when cfg.sim.heteroOn is true, the local rate
%         alpha is multiplied by S.growthField(lin), a smooth static field
%         (arm B). This makes some regions grow faster, bulging into a few
%         large lobes.
%
%   With both off, this function is byte-for-byte Part 2's STEP_PROLIFERATE.

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

displace = isfield(cfg.sim, 'displace') && cfg.sim.displace;
pushRate = 1.0;
if isfield(cfg.sim, 'pushRate'), pushRate = cfg.sim.pushRate; end
heteroOn = isfield(cfg.sim, 'heteroOn') && cfg.sim.heteroOn;

% Occupied orthogonal-neighbour count (0..4) at every occupied site.
kern = [0 1 0; 1 0 1; 0 1 0];
ncount = conv2(double(occ), kern, 'same');

% Frontier = occupied cells with at least one empty neighbour.
empty = ~occ;
if conn == 4
    se = strel('disk', 1, 0);
else
    se = strel('square', 3);
end
frontierMask = occ & imdilate(empty, se);

if displace
    % Second shell: surrounded cells adjacent to the frontier can also divide,
    % via push. (Cells deeper than one shell cannot reach empty space in a
    % single hop, so they are still excluded -- this stays O(perimeter).)
    layer2Mask   = occ & ~frontierMask & imdilate(frontierMask, se);
    eligibleMask = frontierMask | layer2Mask;
else
    eligibleMask = frontierMask;
end

eligible = find(eligibleMask);

if isempty(eligible)
    S.occ = occ; S.age = age; S.clone = clone;
    return
end

order = randperm(numel(eligible));

for k = 1:numel(order)
    lin = eligible(order(k));
    if ~occ(lin), continue; end
    if age(lin) < T_cell, continue; end

    c = ncount(lin);                       % occupied orthogonal neighbours

    a = alphaByClone(clone(lin));
    p = a * (1 - beta * (c / 4));
    if p <= 0, continue; end

    % Growth-rate field modulation (arm B): local alpha *= g(lin).
    if heteroOn && isfield(S, 'growthField') && ~isempty(S.growthField)
        p = p * S.growthField(lin);
    end

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
        elseif displace && rand() < pushRate
            [occ, age, clone, did] = mech.push_displace(lin, occ, age, clone, N, conn, diagW);
            if did
                age(lin) = 0;
                S.totalDivisions = S.totalDivisions + 1;
            end
        end
    end
end

S.occ = occ;
S.age = age;
S.clone = clone;
end
