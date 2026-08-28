function T = fetch_metadata(cfg, verbose)
%FETCH_METADATA Pull ISIC Archive metadata for the cohort collection.
%
%   T = FETCH_METADATA(cfg) returns a table with one row per image in
%   cfg.data.collectionId, carrying the hierarchical diagnosis fields, the
%   acquisition fields and the licensing fields.
%
%   WHY THIS EXISTS -- the ISIC challenge archives ship segmentation masks
%   but no diagnosis labels, so the three-tier grade (nevus / melanoma in
%   situ / invasive melanoma) has to be joined back from the Archive by
%   isic_id. There is no mask endpoint on the API and no diagnosis in the
%   zip, so both halves are mandatory.
%
%   The Archive moved to a hierarchical taxonomy in 2024. The legacy flat
%   fields benign_malignant, mel_class and mel_type are GONE from the live
%   API; diagnosis_1..diagnosis_5 carry that information now, and
%   diagnosis_3 is the level holding the in-situ / invasive split.
%
%   Public collections need no authentication. Paging is cursor-based: the
%   response's 'next' field is a complete URL and is followed verbatim.
%
%   See also IO.BUILD_MANIFEST.

if nargin < 2
    verbose = true;
end

opts = weboptions('ContentType', 'json', 'Timeout', 60, ...
                  'UserAgent', 'MATLAB-Part1-morphology');

url = sprintf('%simages/search/?collections=%d&limit=%d', ...
              cfg.data.apiBase, cfg.data.collectionId, cfg.data.pageLimit);

rows     = {};
pageNum  = 0;
expected = NaN;

while ~isempty(url)
    resp = webread(url, opts);

    if pageNum == 0 && isfield(resp, 'count')
        expected = resp.count;
    end

    results = resp.results;
    % MATLAB returns a struct array only when every record has identical
    % fields. Clinical metadata is sparse, so a cell array is the norm.
    if isstruct(results)
        results = num2cell(results);
    end

    for k = 1:numel(results)
        rows{end+1} = extract_row(results{k}); %#ok<AGROW>
    end

    pageNum = pageNum + 1;
    if verbose
        fprintf('  page %2d  |  %4d / %d records\n', pageNum, numel(rows), expected);
    end

    url = next_url(resp);
end

if isempty(rows)
    T = table();
    return
end

T = struct2table([rows{:}]);

if ~isnan(expected) && height(T) ~= expected
    warning('io:fetch_metadata:countMismatch', ...
        'API reported %d records but %d were retrieved.', expected, height(T));
end
end

% -------------------------------------------------------------------------
function url = next_url(resp)
%NEXT_URL Extract the cursor URL for the following page, or '' when done.
url = '';
if ~isfield(resp, 'next')
    return
end
nxt = resp.next;
if ischar(nxt) || (isstring(nxt) && strlength(nxt) > 0)
    url = char(nxt);
end
end

% -------------------------------------------------------------------------
function r = extract_row(s)
%EXTRACT_ROW Flatten one API record into a table row.

meta     = gf(s, 'metadata', struct());
acq      = gf(meta, 'acquisition', struct());
clin     = gf(meta, 'clinical', struct());

r.isic_id               = str(gf(s, 'isic_id', ''));
r.copyright_license     = str(gf(s, 'copyright_license', ''));
r.attribution           = str(gf(s, 'attribution', ''));

r.image_type            = str(gf(acq, 'image_type', ''));
r.pixels_x              = num(gf(acq, 'pixels_x', NaN));
r.pixels_y              = num(gf(acq, 'pixels_y', NaN));

r.diagnosis_1           = str(gf(clin, 'diagnosis_1', ''));
r.diagnosis_2           = str(gf(clin, 'diagnosis_2', ''));
r.diagnosis_3           = str(gf(clin, 'diagnosis_3', ''));
r.diagnosis_4           = str(gf(clin, 'diagnosis_4', ''));
r.diagnosis_5           = str(gf(clin, 'diagnosis_5', ''));
r.diagnosis_confirm_type = str(gf(clin, 'diagnosis_confirm_type', ''));

r.melanocytic           = logical_or_missing(gf(clin, 'melanocytic', []));
r.mel_thick_mm          = num(gf(clin, 'mel_thick_mm', NaN));
r.mel_ulcer             = logical_or_missing(gf(clin, 'mel_ulcer', []));
r.age_approx            = num(gf(clin, 'age_approx', NaN));
r.sex                   = str(gf(clin, 'sex', ''));
r.anatom_site_1         = str(gf(clin, 'anatom_site_1', ''));
end

% -------------------------------------------------------------------------
function v = gf(s, name, default)
%GF Field access with a default, tolerant of non-struct input.
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function out = str(v)
if isempty(v)
    out = "";
elseif isstring(v) || ischar(v)
    out = string(v);
else
    out = string(v);
end
end

function out = num(v)
if isempty(v)
    out = NaN;
elseif isnumeric(v)
    out = double(v);
else
    % Breslow thickness arrives as a quoted string, e.g. "0.85".
    out = str2double(string(v));
end
end

function out = logical_or_missing(v)
if isempty(v)
    out = NaN;
elseif islogical(v) || isnumeric(v)
    out = double(v);
else
    out = double(strcmpi(string(v), "true"));
end
end
