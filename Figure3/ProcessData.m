%% Load
d = load('PG_3_ex.mat');
ck = d.ck;
 
%% Replace trailing runs of zeros (from an early "break") with a large sentinel value
% If a trial exited its loop early (e.g. via a try/catch break), the
% remaining preallocated entries stay at 0. Detect the first run of 10
% consecutive zeros and mark everything from that point onward as
% diverged (sentinel value 200), rather than treating it as a valid
% (and misleadingly good) cost value.
for trial = 1:size(ck,1)
    row = ck(trial,:);
    for k = 1:numel(row)-10
        if all(row(k:k+9) == 0)
            ck(trial, k:end) = 200;
            break
        end
    end
end
 
%% Replace divergent / negative values with a sentinel value
% Negative sub-optimality gaps indicate numerical issues (should not
% occur analytically), so flag them as diverged with a large value.
ck(ck < 0) = 100;
 
%% Compute mean and standard deviation across trials
ck_mean = mean(ck, 1, 'omitnan');
ck_std  = std(ck, 0, 1, 'omitnan');
save('processed_data.mat', 'ck_mean', 'ck_std');
disp('Done.');
 