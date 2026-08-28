function M = build_manifest(cfg, verbose)
%BUILD_MANIFEST Join images, expert masks and three-tier grades by isic_id.
%
%   M = BUILD_MANIFEST(cfg) returns one row per analysable lesion, with the
%   image path, the expert mask path, the severity grade and the metadata
%   needed for the sensitivity analyses.
%
%   The join is the load-bearing step of Part 1: the challenge archive
%   supplies pixels and masks but no diagnosis, while the Archive API
%   supplies diagnosis but no masks. Only lesions present in BOTH, and whose
%   diagnosis_3 falls in the three-tier scheme, enter the cohort. Melanomas
%   labelled only 'Melanoma, NOS' are excluded on purpose: they carry no
%   in-situ / invasive refinement and would blur the very contrast the study
%   is built on.
%
%   Image and mask dimensions are asserted to match. They come from the same
%   challenge release so a mismatch should never occur, but a silent
%   misalignment would corrupt every descriptor downstream, so it is checked
%   rather than assumed.
%
%   See also IO.FETCH_METADATA, RUN_ANALYSIS.

if nargin < 2
    verbose = true;
end

% ---- Metadata ----------------------------------------------------------
if isfile(cfg.paths.metadata)
    T = readtable(cfg.paths.metadata, 'TextType', 'string');
else
    if verbose
        fprintf('No cached metadata; fetching from the ISIC API...\n');
    end
    T = io.fetch_metadata(cfg, verbose);
    writetable(T, cfg.paths.metadata);
end

% ---- Map diagnosis_3 onto the three severity grades -------------------
gradeMap  = cfg.data.gradeMap;
grade     = strings(height(T), 1);
for k = 1:size(gradeMap, 1)
    grade(T.diagnosis_3 == gradeMap{k, 1}) = gradeMap{k, 2};
end

keep = grade ~= "";
if verbose
    fprintf('Metadata records: %d, of which in the three-tier scheme: %d\n', ...
        height(T), nnz(keep));
end

T     = T(keep, :);
grade = grade(keep);

% ---- Locate image and mask files --------------------------------------
imgDir  = cfg.paths.images;
maskDir = cfg.paths.masks;

assert(isfolder(imgDir),  'Image folder not found: %s', imgDir);
assert(isfolder(maskDir), 'Mask folder not found: %s',  maskDir);

imgPath  = fullfile(imgDir,  T.isic_id + ".jpg");
maskPath = fullfile(maskDir, T.isic_id + cfg.data.maskSuffix);

hasImage = isfile(imgPath);
hasMask  = isfile(maskPath);
present  = hasImage & hasMask;

if verbose
    fprintf('Missing images: %d, missing masks: %d\n', ...
        nnz(~hasImage), nnz(~hasMask));
end

T        = T(present, :);
grade    = grade(present);
imgPath  = imgPath(present);
maskPath = maskPath(present);

% ---- Assemble ----------------------------------------------------------
M = table();
M.isic_id        = T.isic_id;
M.imagePath      = imgPath;
M.maskPath       = maskPath;
M.grade          = categorical(grade, cfg.data.gradeOrder, 'Ordinal', true);
M.diagnosis_3    = T.diagnosis_3;
M.histoConfirmed = T.diagnosis_confirm_type == "histopathology";
M.pixels_x       = T.pixels_x;
M.pixels_y       = T.pixels_y;
M.license        = T.copyright_license;
M.attribution    = T.attribution;
M.melThickMm     = T.mel_thick_mm;

% ---- Verify image and mask geometry agree -----------------------------
dimsOk = true(height(M), 1);
for k = 1:height(M)
    try
        iInfo = imfinfo(char(M.imagePath(k)));
        mInfo = imfinfo(char(M.maskPath(k)));
        dimsOk(k) = iInfo(1).Height == mInfo(1).Height && ...
                    iInfo(1).Width  == mInfo(1).Width;
    catch
        dimsOk(k) = false;
    end
end

if any(~dimsOk)
    warning('io:build_manifest:dimMismatch', ...
        '%d lesions dropped for image/mask dimension mismatch or unreadable files.', ...
        nnz(~dimsOk));
end
M = M(dimsOk, :);

if verbose
    fprintf('\nFinal cohort: %d lesions\n', height(M));
    summary(M.grade);
    fprintf('Histopathology-confirmed: %d of %d\n', ...
        nnz(M.histoConfirmed), height(M));
end

writetable(M, cfg.paths.manifest);
end
