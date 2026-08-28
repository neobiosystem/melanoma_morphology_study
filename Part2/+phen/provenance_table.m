function T = provenance_table()
%PROVENANCE_TABLE Parameter-provenance table for the report.
%
%   T = PROVENANCE_TABLE() returns the P.provenance struct as a two-column
%   table (parameter, source). The benign/melanoma numeric values are joined
%   in so the table reads as "parameter | benign | melanoma | provenance".

P = phen.presets();
names = fieldnames(P.provenance);
rows  = cell(numel(names), 1);
for k = 1:numel(names)
    n = names{k};
    rows{k, 1} = n;
    rows{k, 2} = P.provenance.(n);
end
T = cell2table(rows, 'VariableNames', {'parameter', 'provenance'});
end
