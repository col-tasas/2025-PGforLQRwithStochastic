%% Load and process NPG data
d1 = load('dataNPG_noise1e-5_non_merged.mat');
d2 = load('dataNPG_noise1e-3_non_merged.mat');
%d3 = load('dataNPG_noise1e-1.mat');
d4 = load('dataNPG_noise1e-5_ada_merged.mat');
d5 = load('dataNPG_noise1e-3_ada_merged.mat');
%d6 = load('dataNPG_noise1e-1ada.mat');

m1 = d1.m;
m2 = d2.m;
m4 = d4.m;
m5 = d5.m;

% Ratio: adaptive / fixed * 100 (percentage)
ratio1 = m4 ./ m1 * 100;  % noise 1e-5
ratio2 = m5 ./ m2 * 100;  % noise 1e-3

numPoints = 200;   % approximate number of points to plot per curve (adjust as needed)

% Log-spaced sample indices so the plot has good resolution at both
% early iterations (fast initial transient) and late iterations
% (long-run behavior), without plotting every single data point
idx1 = logIdx(numel(m1), numPoints);
idx2 = logIdx(numel(m2), numPoints);
idx5 = logIdx(numel(m5), numPoints);
idx4 = logIdx(numel(m4), numPoints);

%% Plot
figure; hold on;

% noise 1e-5: adaptive/fixed ratio (m1 and m4 share the same iteration
% count, so idx1 indexes both)
plot(log10(idx1), ratio1(idx1), ...
    'LineStyle','-', 'Color','blue', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-5}I_5$');

% noise 1e-3: adaptive/fixed ratio (m2 and m5 share the same iteration
% count, so idx2 indexes both)
plot(log10(idx2), ratio2(idx2), ...
    'LineStyle','-', 'Color','red', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-3}I_5$');

% reference line at 100% (adaptive == fixed)
yline(100, 'LineStyle','--', 'Color','black', 'LineWidth',1, ...
    'DisplayName', '$100\%$ baseline');

% Manually set x-axis ticks to display original (non-log) values, since
% the x-axis data itself is plotted on a log10 scale
ax = gca;
ax.XAxis.TickValues = log10([1, 10, 100, 1000, 10000, 100000, 1000000]);
ax.XAxis.TickLabels = {'1','10','100','10^3','10^4','10^5','10^6'};
xlim([log10(1000),log10(275000)])

grid on;
box on;
ylabel({'Suboptimality ratio:' '$\frac{\mathrm{adaptive}}{\mathrm{fixed}}$ (\%)'},'interpreter','latex','FontSize',12);
xlabel('$\mathrm{iteration}~i$','interpreter','latex','FontSize',12);
legend('interpreter','latex','FontSize',12,'Location','best');
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