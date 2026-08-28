function [nbr, ok] = sample_empty_neighbor(lin, occ, N, conn, diagWeight, C, chi)
%SAMPLE_EMPTY_NEIGHBOR Pick one empty neighbour of a site, isotropically.
%
%   [nbr, ok] = SAMPLE_EMPTY_NEIGHBOR(lin, occ, N, conn, diagWeight, C, chi)
%   samples an empty neighbour of the occupied site LIN. Weights come from
%   CA.NEIGHBORS, optionally multiplied by exp(chi * (C_neighbour - C_self))
%   for a chemotactic/haptotactic bias up the nutrient gradient.
%
%   Returns nbr = 0 and ok = false when the site has no empty neighbour.

[nidx, w] = ca.neighbors(lin, N, conn, diagWeight);

empty = ~occ(nidx);
nidx = nidx(empty);
w    = w(empty);

if isempty(nidx)
    nbr = 0; ok = false;
    return
end

if nargin >= 7 && ~isempty(C) && chi ~= 0
    w = w .* exp(chi * (C(nidx) - C(lin)));
end

w = w / sum(w);
u = rand();
cs = cumsum(w);
k = find(u <= cs, 1, 'first');

nbr = nidx(k);
ok  = true;
end
