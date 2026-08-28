function S = init_state(cfg)
%INIT_STATE Build the initial state of a CA simulation.
%
%   S = INIT_STATE(cfg) returns a struct S holding the occupancy lattice
%   (logical N x N), per-cell age and clone id (both 0 off-lattice), and the
%   nutrient-field placeholders. The seed is a disk of initRadius sites at
%   the frame centre; a fraction subclone.frac of its cells are assigned to
%   the advantaged clone (mechanism arm D).

N = cfg.sim.N;
occ = false(N, N);

c = round(N / 2);
[XX, YY] = meshgrid(1:N, 1:N);
occ = (XX - c).^2 + (YY - c).^2 <= cfg.sim.initRadius.^2;

age   = zeros(N, N);
clone = int32(occ);                 % clone 1 where occupied

if cfg.sim.subclone.frac > 0
    idx  = find(occ);
    k    = max(1, round(cfg.sim.subclone.frac * numel(idx)));
    pick = idx(randperm(numel(idx), k));
    clone(pick) = 2;
end

S = struct();
S.cfg = cfg;
S.N   = N;
S.occ = occ;
S.age = age;
S.clone = clone;

S.C     = [];                       % coarse nutrient field (nc x nc)
S.Cfine = [];                       % fine interpolation (N x N)
S.rho   = [];                       % coarse occupancy fraction

S.step   = 0;
S.nCells = nnz(occ);
S.initCells = nnz(occ);
S.totalDivisions = 0;
S.totalDeaths    = 0;

% Per-clone proliferation rate; clone 1 is the parent clone, clone 2 the
% advantaged subclone (mechanism arm D).
S.alphaByClone = [cfg.sim.alpha; cfg.sim.alpha * cfg.sim.subclone.alphaMult];
end
