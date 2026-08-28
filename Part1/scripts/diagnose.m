function diagnose()
%DIAGNOSE Print raw numbers behind the failing validation cases.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root,'config')); addpath(fullfile(root,'tests'));
cfg = default_config();

show = @(name, o) fprintf(['%-22s Df=%7.4f  R2=%7.4f  n=%3d  span=%.2f  ' ...
    'eps=[%7.2f %7.2f]  flag=%s\n'], name, o.Df, o.R2, o.numPoints, ...
    o.spanDecades, o.epsRange(1), o.epsRange(2), o.flag);

fprintf('=== box counting ===\n');

P = shapes.diagonal_line(800);
o = morph.fractal_boxcount(P, 800, cfg);
show('diagonal line', o);
fprintf('   npts=%d  counts(1:5)=%s\n', size(P,1), mat2str(o.countsAll(1:5)'));

m = shapes.square(600); P = morph.boundary_pixels(m);
d = morph.equivalent_diameter(m);
o = morph.fractal_boxcount(P, d, cfg);
show('square boundary', o);
fprintf('   dEq=%.1f  npts=%d  counts(1:5)=%s\n', d, size(P,1), mat2str(o.countsAll(1:5)'));

P = shapes.filled_disk_points(400);
o = morph.fractal_boxcount(P, 800, cfg);
show('filled disk', o);
fprintf('   npts=%d  counts=%s\n', size(P,1), mat2str(round(o.countsAll(1:8))'));
fprintf('   local slopes=%s\n', mat2str(round(o.localSlopes(~isnan(o.localSlopes))',3)));

[P,d] = shapes.koch_points(5, 2048);
o = morph.fractal_boxcount(P, d, cfg);
show('koch order5', o);
fprintf('   dEq=%.1f  npts=%d\n', d, size(P,1));
fprintf('   local slopes=%s\n', mat2str(round(o.localSlopes(~isnan(o.localSlopes))',3)));

fprintf('\n=== compactness of analytic circle ===\n');
try
    C = shapes.circle_contour(400);
    fprintf('   size(C)=%s  first=%s last=%s\n', mat2str(size(C)), ...
        mat2str(round(C(1,:),3)), mat2str(round(C(end,:),3)));
    fprintf('   perim=%.3f (true %.3f)  area=%.1f (true %.1f)\n', ...
        morph.contour_perimeter(C), 2*pi*400, morph.contour_area(C), pi*400^2);
    fprintf('   compactness=%.6f\n', morph.compactness(C));
catch ME
    fprintf('   ERROR: %s\n   %s\n', ME.identifier, ME.message);
    for k=1:min(3,numel(ME.stack))
        fprintf('     at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
    end
end

fprintf('\n=== half disk asymmetry ===\n');
m = shapes.disk(300);
m(1:round(size(m,1)/2), :) = false;
a = morph.asymmetry_index(m);
fprintf('   major=%.4f minor=%.4f mean=%.4f central=%.4f elong=%.4f\n', ...
    a.major, a.minor, a.mean, a.central, a.elongation);
end
