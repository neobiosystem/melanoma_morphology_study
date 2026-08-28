classdef shapes
%SHAPES Synthetic test geometry with analytically known properties.
%
%   Generators used by the Part 1 validation suite. Every shape here has a
%   closed-form fractal dimension, perimeter or compactness, so a failing
%   test points at the implementation rather than at the data.

methods (Static)

    function mask = disk(radius, pad)
        %DISK Filled circular mask of the given radius.
        if nargin < 2, pad = 20; end
        n = 2 * (radius + pad) + 1;
        c = radius + pad + 1;
        [X, Y] = meshgrid(1:n, 1:n);
        mask = hypot(X - c, Y - c) <= radius;
    end

    function mask = square(side, pad)
        %SQUARE Filled square mask of the given side length.
        if nargin < 2, pad = 20; end
        n = side + 2 * pad;
        mask = false(n, n);
        mask(pad + (1:side), pad + (1:side)) = true;
    end

    function mask = ellipse(semiMajor, semiMinor, pad)
        %ELLIPSE Filled axis-aligned ellipse mask.
        if nargin < 3, pad = 20; end
        h = 2 * (semiMinor + pad) + 1;
        w = 2 * (semiMajor + pad) + 1;
        cy = semiMinor + pad + 1;
        cx = semiMajor + pad + 1;
        [X, Y] = meshgrid(1:w, 1:h);
        mask = ((X - cx) / semiMajor).^2 + ((Y - cy) / semiMinor).^2 <= 1;
    end

    function P = circle_points(radius, numPts)
        %CIRCLE_POINTS Analytic circle boundary, as [row col] coordinates.
        if nargin < 2, numPts = 4000; end
        t = linspace(0, 2*pi, numPts + 1)';
        t(end) = [];
        P = [radius * sin(t), radius * cos(t)];
    end

    function C = circle_contour(radius, numPts)
        %CIRCLE_CONTOUR Closed analytic circle contour.
        if nargin < 2, numPts = 4000; end
        P = shapes.circle_points(radius, numPts);
        C = [P; P(1, :)];
    end

    function C = square_contour(side, ptsPerEdge)
        %SQUARE_CONTOUR Closed analytic square contour, densely sampled.
        if nargin < 2, ptsPerEdge = 1000; end
        t = linspace(0, 1, ptsPerEdge + 1)';
        t(end) = [];
        edge1 = [zeros(size(t)),        t * side];
        edge2 = [t * side,              side * ones(size(t))];
        edge3 = [side * ones(size(t)),  side - t * side];
        edge4 = [side - t * side,       zeros(size(t))];
        C = [edge1; edge2; edge3; edge4];
        C = [C; C(1, :)];
    end

    function V = koch_vertices(order)
        %KOCH_VERTICES Vertices of the Koch snowflake at the given order.
        %   The snowflake boundary has fractal dimension log(4)/log(3).
        V = [0, 0; 1, 0; 0.5, sqrt(3)/2; 0, 0];

        for iter = 1:order
            nSeg = size(V, 1) - 1;
            newV = zeros(4 * nSeg + 1, 2);

            for k = 1:nSeg
                p1 = V(k, :);
                p2 = V(k + 1, :);
                d  = (p2 - p1) / 3;
                a  = p1 + d;
                b  = p1 + 2 * d;

                % Rotate the middle third by -60 deg so the bump points
                % outward for a counter-clockwise triangle.
                ang  = -pi / 3;
                R    = [cos(ang), -sin(ang); sin(ang), cos(ang)];
                apex = a + (R * d')';

                newV(4*(k-1) + (1:4), :) = [p1; a; apex; b];
            end

            newV(end, :) = V(end, :);
            V = newV;
        end
    end

    function [P, dEq] = koch_points(order, sizePx)
        %KOCH_POINTS Rasterised Koch snowflake boundary and its D_eq.
        %   Returns the digital boundary point set together with the
        %   equivalent diameter of the enclosed region, which the box
        %   counter needs to set its scale bounds.
        V = shapes.koch_vertices(order);
        V = V * sizePx;
        V = V - min(V, [], 1) + 1;

        P = shapes.rasterize_polyline(V, 0.5);

        area = polyarea(V(1:end-1, 1), V(1:end-1, 2));
        dEq  = 2 * sqrt(area / pi);
    end

    function P = rasterize_polyline(V, step)
        %RASTERIZE_POLYLINE Sample a polyline onto unique integer pixels.
        if nargin < 2, step = 0.5; end

        nSeg    = size(V, 1) - 1;
        segLen  = hypot(diff(V(:,1)), diff(V(:,2)));
        nPer    = max(2, ceil(segLen / step));
        total   = sum(nPer);

        P = zeros(total, 2);
        pos = 1;
        for k = 1:nSeg
            n = nPer(k);
            t = linspace(0, 1, n)';
            P(pos:pos+n-1, :) = V(k, :) + t .* (V(k+1, :) - V(k, :));
            pos = pos + n;
        end

        P = unique(round(P), 'rows');
    end

    function P = diagonal_line(len)
        %DIAGONAL_LINE Digital diagonal line; dimension 1.
        t = (0:0.5:len)';
        P = unique(round([t, t]), 'rows');
    end

    function P = filled_disk_points(radius)
        %FILLED_DISK_POINTS Every pixel of a disk; dimension 2.
        mask = shapes.disk(radius, 5);
        [r, c] = find(mask);
        P = [r, c];
    end

end
end
