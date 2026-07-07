%% Load data
% Load results from the four PG experiments (see policy_gradient_lqr.m)
d1 = load('dataPG_noisefree.mat');
d2 = load('dataPG_noiselow_stephigh.mat');
d3 = load('dataPG_noisehigh_stephigh.mat');
d4 = load('dataPG_noisehigh_steplow.mat');
 
step = 500;        % subsampling interval (unused now that log-spaced indices are used, kept for reference)
numPoints = 200;    % approximate number of points to plot per curve (adjust as needed)
 
% Log-spaced sample indices so the plot has good resolution at both
% early iterations (fast initial transient) and late iterations
% (long-run behavior), without plotting every single data point
idx1 = logIdx(numel(d1.CKPGV), numPoints);
idx2 = logIdx(numel(d2.ck_mean), numPoints);
idx4 = logIdx(numel(d4.ck_mean), numPoints);
 
%% Plot
figure; hold on;
 
% --- Noise-free PG ---
plot(log10(idx1), log10(d1.CKPGV(idx1)), ...
    'LineStyle', ':', 'Color', 'black', 'LineWidth', 1.5, ...
    'DisplayName', '$\sigma^2=0,~\eta=0.0005$');
 
% --- Low noise variance, large step size ---
m2 = d2.ck_mean(idx2); s2 = d2.ck_std(idx2);
% shaded region shows mean +/- 1 std across Monte Carlo runs
patch([log10(idx2), fliplr(log10(idx2))], ...
      [log10(max(m2+s2, 1e-10)), fliplr(log10(max(m2-s2, 1e-10)))], ...
    'blue', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(log10(idx2), log10(m2), ...
    'LineStyle', '--', 'Color', 'blue', 'LineWidth', 1.5, ...
    'DisplayName', ['$\sigma^2=', num2str(d2.noisevar_low), ',~\eta=', num2str(d2.step_high), '$']);
 
% --- High noise variance, large step size ---
% Only the first 16*step samples are used here, since this run tends to
% become unstable/diverge at later iterations
idx3 = 1:step:min(16*step, numel(d3.ck_mean));
m3 = d3.ck_mean(idx3); s3 = d3.ck_std(idx3);
patch([log10(idx3), fliplr(log10(idx3))], ...
      [log10(max(m3+s3, 1e-10)), fliplr(log10(max(m3-s3, 1e-10)))], ...
    'green', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(log10(idx3), log10(m3), ...
    'LineStyle', '-.', 'Color', 'green', 'LineWidth', 1.5, ...
    'DisplayName', ['$\sigma^2=', num2str(d3.noisevar_high), ',~\eta=', num2str(d3.step_high), '$']);
 
% --- High noise variance, small step size ---
m4 = d4.ck_mean(idx4); s4 = d4.ck_std(idx4);
patch([log10(idx4), fliplr(log10(idx4))], ...
      [log10(max(m4+s4, 1e-10)), fliplr(log10(max(m4-s4, 1e-10)))], ...
    'red', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(log10(idx4), log10(m4), ...
    'LineStyle', '-', 'Color', 'red', 'LineWidth', 1.5, ...
    'DisplayName', ['$\sigma^2=', num2str(d4.noisevar_high), ',~\eta=', num2str(d4.step_low), '$']);
 
% Manually set axis ticks to display original (non-log) values, since the
% data itself is plotted on a log10 scale
ax = gca;
ax.XAxis.TickValues = log10([1, 10, 100, 1000, 10000, 100000, 1000000]);
ax.XAxis.TickLabels = {'1', '10', '100', '10^3', '10^4', '10^5', '10^6'};
ax.YAxis.TickValues = log10([0.001, 0.01, 0.1, 1, 10]);
ax.YAxis.TickLabels = {'10^{-3}', '10^{-2}', '10^{-1}', '10^0', '10^1', '10^2'};
ylim([log10(0.001), log10(10)]);
xlim([log10(1), log10(500000)]);
 
grid on;
box on;
ylabel('$\frac{C({K}_i)-C(K^*)}{C(K^*)}$', 'interpreter', 'latex', 'FontSize', 12);
xlabel('$\mathrm{iteration}~i$', 'interpreter', 'latex', 'FontSize', 12);
legend('interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
hold off;
 
function idx = logIdx(N, numPoints)
% logIdx  Generate log-spaced, unique, sorted integer indices in [1, N]
%   idx = logIdx(N, numPoints) returns approximately numPoints indices
%   spaced evenly on a log scale between 1 and N, useful for subsampling
%   long time series for plotting without losing early-time detail.
    idx = round(logspace(0, log10(N), numPoints));
    idx = unique(idx);
    idx(idx < 1) = 1;
    idx(idx > N) = N;
end