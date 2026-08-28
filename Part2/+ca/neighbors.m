function [idx, w] = neighbors(lin, N, conn, diagWeight)
%NEIGHBORS Valid in-bounds neighbour indices of one lattice site.
%
%   [idx, w] = NEIGHBORS(lin, N, conn, diagWeight) returns the linear
%   indices of the von Neumann (conn = 4) or Moore (conn = 8) neighbours of
%   the site LIN, together with their isotropy weights W.
%
%   Moore diagonals are down-weighted by diagWeight = 1/sqrt(2) so that a
%   growth or migration step is isotropic: a diagonal hop is sqrt(2) further,
%   so it must be chosen sqrt(2) times less often. Without this correction a
%   square lattice grows a spurious 4-fold (diamond/square) artefact, which
%   would show up as a fake compactness difference -- the exact failure
%   TestIsotropy exists to catch.

if nargin < 3 || isempty(conn), conn = 4; end
if nargin < 4 || isempty(diagWeight), diagWeight = 1/sqrt(2); end

[r, c] = ind2sub([N N], lin);

if conn == 4
    dr = [-1  1  0  0];
    dc = [ 0  0 -1  1];
    w  = [1 1 1 1];
else
    dr = [-1  1  0  0 -1 -1  1  1];
    dc = [ 0  0 -1  1 -1  1 -1  1];
    w  = [1 1 1 1 diagWeight diagWeight diagWeight diagWeight];
end

rr = r + dr;
cc = c + dc;
ok = rr >= 1 & rr <= N & cc >= 1 & cc <= N;

idx = sub2ind([N N], rr(ok), cc(ok));
w   = w(ok);
end
