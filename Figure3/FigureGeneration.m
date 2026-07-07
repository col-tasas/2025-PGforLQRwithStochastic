%% Load and process NPG data
d1 = load('PG_5.mat');
d2 = load('PG_3.mat');
d3 = load('PG_VR5.mat');
d4 = load('PG_VR3.mat');
d5 = load('processed_datatest.mat');

% Mean and std
m1 = mean(d1.ck, 1, 'omitnan');  s1 = std(d1.ck, 0, 1, 'omitnan');
m2 = mean(d2.ck, 1, 'omitnan');  s2 = std(d2.ck, 0, 1, 'omitnan');
m3 = mean(d3.ck, 1, 'omitnan');  s3 = std(d3.ck, 0, 1, 'omitnan');
m4 = mean(d4.ck, 1, 'omitnan');  s4 = std(d4.ck, 0, 1, 'omitnan');
m5 = abs(d5.ck_mean);             s5 = d5.ck_std;

numPoints = 200;   % approximate number of points to plot per curve (adjust as needed)

% Log-spaced sample indices so the plot has good resolution at both
% early iterations (fast initial transient) and late iterations
% (long-run behavior), without plotting every single data point
idx1 = logIdx(numel(m1), numPoints);
idx2 = logIdx(numel(m2), numPoints);
idx3 = logIdx(numel(m3), numPoints);
idx4 = logIdx(numel(m4), numPoints);

%% Plot
figure; hold on;

% --- Sigma_w = 1e-5*I, eta = 0.1 (no variance reduction) ---
patch([log10(idx1), fliplr(log10(idx1))], ...
      [log10(max(m1(idx1)+s1(idx1), 1e-10)), fliplr(log10(max(m1(idx1)-s1(idx1), 1e-10)))], ...
    'blue', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(log10(idx1), log10(m1(idx1)), ...
    'LineStyle','--', 'Color','blue', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-5}I_5,~\eta=0.1$');

% --- processed_datatest.mat curve (only first 100 points, step 10) ---
idx5 = 1:10:100;
patch([log10(idx5), fliplr(log10(idx5))], ...
      [log10(max(m5(idx5)+s5(idx5), 1e-10)), fliplr(log10(max(m5(idx5)-s5(idx5), 1e-10)))], ...
    'black', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(log10(idx5), log10(m5(idx5)), ...
    'LineStyle','-.', 'Marker','o','Color','black', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-3}I_5,~\eta=0.002$');   % update DisplayName if the underlying settings change

% --- Sigma_w = 1e-3*I, eta = 0.0008 (no variance reduction) ---
patch([log10(idx2), fliplr(log10(idx2))], ...
      [log10(max(m2(idx2)+s2(idx2), 1e-10)), fliplr(log10(max(m2(idx2)-s2(idx2), 1e-10)))], ...
    'green', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(log10(idx2), log10(m2(idx2)), ...
    'LineStyle',':', 'Color','green', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-3}I_5,~\eta=0.0008$');

% --- Sigma_w = 1e-5*I, eta = 0.1, with variance reduction (VR) ---
patch([log10(idx3), fliplr(log10(idx3))], ...
      [log10(max(m3(idx3)+s3(idx3), 1e-10)), fliplr(log10(max(m3(idx3)-s3(idx3), 1e-10)))], ...
    'red', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(log10(idx3), log10(m3(idx3)), ...
    'LineStyle','-', 'Color','red', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-5}I_5,~\eta=0.1~\mathrm{VR}$');

% --- Sigma_w = 1e-3*I, eta = 0.0008, with variance reduction (VR) ---
patch([log10(idx4), fliplr(log10(idx4))], ...
      [log10(max(m4(idx4)+s4(idx4), 1e-10)), fliplr(log10(max(m4(idx4)-s4(idx4), 1e-10)))], ...
    'magenta', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(log10(idx4), log10(m4(idx4)), ...
    'LineStyle',':', 'Color','magenta', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-3}I_5,~\eta=0.0008~\mathrm{VR}$');

% Manually set axis ticks to display original (non-log) values, since the
% data itself is plotted on a log10 scale
ax = gca;
ax.XAxis.TickValues = log10([0, 1, 10, 100, 1000, 10000, 100000, 1000000]);
ax.XAxis.TickLabels = {'1','10','100','10^3','10^4','10^5','10^6','10^7'};
ax.YAxis.TickValues = log10([0.001, 0.01, 0.1, 1, 10, 100]);
ax.YAxis.TickLabels = {'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'};
ylim([log10(0.0005), log10(3)]);
xlim([log10(1), log10(3500000)]);

grid on;
box on;
ylabel('$\frac{C(\hat{K}_i)-C(K^*)}{C(K^*)}$','interpreter','latex','FontSize',12);
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