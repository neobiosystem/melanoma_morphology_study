classdef TestSegmentation < matlab.unittest.TestCase
%TESTSEGMENTATION Validation of the automatic segmentation route.
%
%   The tests that matter here concern vignette handling, because it is the
%   one place where a plausible-looking rule silently returns the wrong
%   object rather than failing loudly.
%
%   Dermatoscopes image through a circular aperture, so frame corners are
%   dark and Otsu labels them 'lesion'. The obvious fix -- drop everything
%   touching the border -- is wrong for a lesion that genuinely runs off the
%   frame: any stray interior speck then becomes the answer and the lesion is
%   discarded. Both directions are tested, since a rule that passes only one
%   is exactly the failure mode being guarded against.

properties
    cfg
end

methods (TestClassSetup)
    function loadConfig(testCase)
        c = default_config();
        % Active-contour refinement is a separate concern; disabling it keeps
        % these tests measuring the component-selection logic alone.
        c.seg.useActiveContour = false;
        testCase.cfg = c;
    end
end

methods (Test)

    function testVignetteCornersAreDiscarded(testCase)
        % Dark corners plus a clearly interior lesion: the corners must go.
        [img, lesion] = TestSegmentation.scene_with_vignette();
        mask = seg.otsu_route(img, testCase.cfg);

        recall = nnz(mask & lesion) / nnz(lesion);
        testCase.verifyGreaterThan(recall, 0.85, ...
            'The interior lesion must be recovered.');

        corners = false(size(lesion));
        corners([1:120, end-119:end], [1:120, end-119:end]) = true;
        testCase.verifyLessThan(nnz(mask & corners) / max(nnz(mask), 1), 0.05, ...
            'Vignette corners must not survive into the final mask.');
    end

    function testBorderTouchingLesionIsNotReplacedBySpeck(testCase)
        % The regression this guards: a lesion running off the left edge plus
        % one small interior speck. Clearing all border components would
        % return the speck, which scores a near-zero Dice against the expert
        % mask while reporting itself as a valid segmentation.
        [img, lesion, speck] = TestSegmentation.scene_with_border_lesion();
        mask = seg.otsu_route(img, testCase.cfg);

        recall = nnz(mask & lesion) / nnz(lesion);
        testCase.verifyGreaterThan(recall, 0.85, ...
            'A lesion touching the frame must be kept, not cleared.');
        testCase.verifyLessThan(nnz(mask & speck) / max(nnz(mask), 1), 0.10, ...
            'The interior speck must not become the selected component.');
    end

    function testSingleComponentIsReturned(testCase)
        [img, ~] = TestSegmentation.scene_with_vignette();
        mask = seg.otsu_route(img, testCase.cfg);
        cc = bwconncomp(mask, 8);
        testCase.verifyEqual(cc.NumObjects, 1, ...
            'The route must commit to exactly one lesion component.');
    end

    function testLesionIsDarkerClassNotLighter(testCase)
        % Polarity check. If the threshold comparison were inverted the route
        % would return the skin and every descriptor would describe the
        % background instead of the lesion.
        [img, lesion] = TestSegmentation.scene_with_vignette();
        mask = seg.otsu_route(img, testCase.cfg);
        testCase.verifyLessThan(nnz(mask), nnz(~lesion), ...
            'The mask must be the dark lesion, not the surrounding skin.');
    end

end

methods (Static)

    function [img, lesion] = scene_with_vignette()
        %SCENE_WITH_VIGNETTE Interior dark lesion plus dark frame corners.
        n = 600;
        [X, Y] = meshgrid(1:n, 1:n);

        lesion = (X - n/2).^2 + (Y - n/2).^2 <= 120^2;

        corners = false(n);
        for cx = [1, n]
            for cy = [1, n]
                corners = corners | ((X - cx).^2 + (Y - cy).^2 <= 100^2);
            end
        end

        img = TestSegmentation.render(n, lesion | corners);
    end

    function [img, lesion, speck] = scene_with_border_lesion()
        %SCENE_WITH_BORDER_LESION Lesion running off the left edge, plus a
        %small interior speck that must not win.
        n = 600;
        [X, Y] = meshgrid(1:n, 1:n);

        lesion = false(n);
        lesion(100:500, 1:250) = true;

        speck = (X - 450).^2 + (Y - 300).^2 <= 22^2;

        img = TestSegmentation.render(n, lesion | speck);
    end

    function img = render(n, dark)
        %RENDER Build an RGB scene: light skin, dark lesion, mild noise.
        % Noise is seeded so a failure is reproducible rather than
        % intermittent.
        rng(7);
        base = uint8(200 * ones(n, n));
        base(dark) = 60;
        base = uint8(double(base) + 4 * randn(n));
        img  = repmat(base, 1, 1, 3);
    end

end
end
