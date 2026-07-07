%% ========================================================================
%  Clean and summarize NPG data
%  For each of the four files below: remove any trial (row) containing
%  NaN or any value > 10, then compute the mean and std across the
%  remaining trials.
%  ========================================================================
clear; clc;

files = {'NPG5Ada.mat', 'NPG3Ada.mat', 'NPG5NonAda.mat', 'NPG3NonAda.mat'};

for iFile = 1:length(files)
    fname = files{iFile};
    d = load(fname);
    ck = d.ck;   % rows = Monte Carlo trials, columns = iterations

    % remove any trial (entire row) that contains NaN or a value > 10
    badRows = any(isnan(ck), 2) | any(ck > 10, 2);
    nRemoved = sum(badRows);
    ck(badRows, :) = [];

    fprintf('%s: %d trials before, %d removed (NaN or >10), %d remaining\n', ...
        fname, nRemoved + size(ck,1), nRemoved, size(ck,1));

    % mean / std across the remaining trials
    m = mean(ck, 1, 'omitnan');
    s = std(ck, 0, 1, 'omitnan');

    % save with a name derived from the input file, e.g. NPG5Ada_processed.mat
    [~, baseName] = fileparts(fname);
    outName = sprintf('%s_processed.mat', baseName);
    save(outName, 'ck', 'm', 's');
    fprintf('  -> saved to %s\n\n', outName);
end

disp('All files processed.');