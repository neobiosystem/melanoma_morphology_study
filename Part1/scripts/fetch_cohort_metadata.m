function T = fetch_cohort_metadata()
%FETCH_COHORT_METADATA Download and summarise the ISIC cohort metadata.
%
%   Pulls every record in the cohort collection from the ISIC Archive API,
%   writes it to data/isic_metadata.csv and prints the grade breakdown.
%   Run this once before building the manifest.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root, 'config'));

cfg = default_config();

fprintf('Fetching metadata for collection %d ...\n', cfg.data.collectionId);
T = io.fetch_metadata(cfg, true);

fprintf('\nRetrieved %d records, %d fields.\n\n', height(T), width(T));

fprintf('--- diagnosis_3 (the three-tier grading field) ---\n');
disp(sortrows(groupcounts(T, 'diagnosis_3'), 'GroupCount', 'descend'));

fprintf('--- image_type ---\n');
disp(groupcounts(T, 'image_type'));

fprintf('--- diagnosis_confirm_type ---\n');
disp(groupcounts(T, 'diagnosis_confirm_type'));

fprintf('--- copyright_license ---\n');
disp(groupcounts(T, 'copyright_license'));

writetable(T, cfg.paths.metadata);
fprintf('Wrote %s\n', cfg.paths.metadata);
end
