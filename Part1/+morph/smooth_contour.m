function C = smooth_contour(C0, windowPx)
%SMOOTH_CONTOUR Periodic moving-average smoothing of a closed contour.
%
%   C = SMOOTH_CONTOUR(C0, windowPx) smooths the closed contour C0 (N-by-2,
%   first vertex repeated last) with a moving average of WINDOWPX vertices
%   and returns a contour that is still explicitly closed.
%
%   Why this exists: perimeter measured on a raster boundary is inflated by
%   staircase effects, and inflated MORE for irregular boundaries than for
%   smooth ones. Left uncorrected that bias would manufacture a spurious
%   between-grade difference in compactness -- melanomas would look less
%   compact purely because their perimeter was over-measured. Callers pass a
%   window that is a fixed fraction of the lesion's equivalent diameter, not
%   a fixed pixel count, so the correction is resolution-independent.
%
%   The contour is triplicated before filtering and the middle copy is
%   returned. A plain movmean would shrink the window at the endpoints and
%   leave a visible seam where the contour closes on itself.
%
%   See also MORPH.CONTOUR_PERIMETER, MORPH.COMPACTNESS.

% Work on the open vertex list; closure is restored at the end.
if size(C0, 1) > 1 && isequal(C0(1, :), C0(end, :))
    C0(end, :) = [];
end

n = size(C0, 1);
if n < 3
    C = [C0; C0(1, :)];
    return
end

windowPx = 2 * floor(windowPx / 2) + 1;   % force odd so the filter is centred
if windowPx <= 1
    C = [C0; C0(1, :)];
    return
end
windowPx = min(windowPx, n);

replicated = [C0; C0; C0];
smoothed   = movmean(replicated, windowPx, 1);
C          = smoothed(n + 1 : 2 * n, :);
C          = [C; C(1, :)];
end
