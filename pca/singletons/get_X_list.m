function [Xlist, Vlist, grpSize] = get_X_list(numComponents)

subjects = [1,3,4,9,17,35,36,37,39,44,48,50,51,52,53,54,55,66,69,71,75,76,78,79,81];

G = 18;
Xlist = cell(1, numel(subjects));
Vlist = cell(1, numel(subjects));
grpSize = zeros(1, G);

for N = 1:numel(subjects)
    Xlist{N} = cell(1, G);
    load(['forwardAndRois-skeri' num2str(subjects(N))]);
    load(['ROI_correlation_subj_' num2str(subjects(N))]);
    for g = 1:G
        grpSize(g) = grpSize(g) + numel(ROIs.ndx{g});
        [Xlist{N}{g}, Vlist{N}{g}] = get_principal_components(fwdMatrix(:, ROIs.ndx{g}), numComponents);
    end
end

end