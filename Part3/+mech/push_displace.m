function [occ, age, clone, did] = push_displace(lin, occ, age, clone, N, conn, diagW)
%PUSH_DISPLACE Local O(1) push: displace one occupied neighbour outward.
%
%   [occ, age, clone, did] = PUSH_DISPLACE(lin, occ, age, clone, N, conn, diagW)
%   is the mechanical resolution of a surrounded cell's division. LIN has no
%   empty neighbour, so instead of dropping the division we pick one of its
%   occupied neighbours NB (isotropically) and move NB into one of NB's own
%   empty neighbours, then place the daughter at NB. Net +1 cell, so
%   nCells = initCells + totalDivisions - totalDeaths still holds.
%
%   This is a single-hop push (the neighbour moves one site), which resolves
%   the "surrounded cell cannot divide" bottleneck one layer in from the
%   frontier at O(1) -- not a full-chain relaxation to the free surface.

did = false;

[nidx, w] = ca.neighbors(lin, N, conn, diagW);
occNb = occ(nidx);
nidx  = nidx(occNb);
w     = w(occNb);
if isempty(nidx), return; end

w = w / sum(w);
u = rand();
cs = cumsum(w);
k = find(u <= cs, 1, 'first');
nb = nidx(k);

[target, ok] = ca.sample_empty_neighbor(nb, occ, N, conn, diagW, [], 0);
if ~ok, return; end

% Move neighbour nb -> target, daughter occupies nb (read nb's identity first).
occ(target)   = true;
age(target)   = age(nb);
clone(target) = clone(nb);
occ(nb)       = true;
age(nb)       = 0;
clone(nb)     = clone(lin);

did = true;
end
