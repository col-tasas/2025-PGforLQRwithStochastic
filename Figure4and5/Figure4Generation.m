%% Load and process NPG data
d1 = load('NPG5NonAda_processed.mat');   % Sigma_w = 1e-5, fixed step size
d2 = load('NPG3NonAda_processed.mat');   % Sigma_w = 1e-3, fixed step size
%d3 = load('dataNPG_noise1e-1.mat');      % Sigma_w = 1e-1, fixed step size (unused)
d4 = load('NPG5Ada_processed.mat');      % Sigma_w = 1e-5, adaptive step size
d5 = load('NPG3Ada_processed.mat');      % Sigma_w = 1e-3, adaptive step size
%d6 = load('dataNPG_noise1e-1ada.mat');   % Sigma_w = 1e-1, adaptive step size (unused)

% mean and std were already computed when the data was processed/cleaned
m1 = d1.m;  s1 = d1.s;
m2 = d2.m;  s2 = d2.s;
%m3 = mean(d3.ck, 1, 'omitnan');  s3 = std(d3.ck, 0, 1, 'omitnan');
m4 = d4.m;  s4 = d4.s;
m5 = d5.m;  s5 = d5.s;
%m6 = mean(d6.ck, 1, 'omitnan');  s6 = std(d6.ck, 0, 1, 'omitnan');

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

% --- Sigma_w = 1e-5, fixed step size ---
patch([log10(idx1), fliplr(log10(idx1))], ...
      [log10(max(m1(idx1)+s1(idx1), 1e-10)), fliplr(log10(max(m1(idx1)-s1(idx1), 1e-10)))], ...
      'blue', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(log10(idx1), log10(m1(idx1)), ...
    'LineStyle',':', 'Color','blue', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-5}~\mathrm{Fixed}~\eta_i=\frac{1}{2+25000\mathrm{Tr}(P_{\hat{K}_0})}$');

% --- Sigma_w = 1e-3, fixed step size ---
patch([log10(idx2), fliplr(log10(idx2))], ...
      [log10(max(m2(idx2)+s2(idx2), 1e-10)), fliplr(log10(max(m2(idx2)-s2(idx2), 1e-10)))], ...
      'red', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(log10(idx2), log10(m2(idx2)), ...
    'LineStyle','-.', 'Color','red', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-3}~\mathrm{Fixed}~\eta_i=\frac{1}{2+25000\mathrm{Tr}(P_{\hat{K}_0})}$');

% --- Sigma_w = 1e-1, fixed step size (currently disabled) ---
% idx3 = 1:step:numel(m3);
% patch([log10(idx3), fliplr(log10(idx3))], ...
%       [log10(max(m3(idx3)+s3(idx3), 1e-10)), fliplr(log10(max(m3(idx3)-s3(idx3), 1e-10)))], ...
%       'red', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
% plot(log10(idx3), log10(m3(idx3)), ...
%     'LineStyle','-', 'Color','red', 'LineWidth',1.5, ...
%     'DisplayName', '$\Sigma_w=10^{-1}~\eta_i=\frac{1}{2+25000\mathrm{Tr}(P_{\hat{K}_0})}$');

% --- Sigma_w = 1e-5, adaptive step size ---
patch([log10(idx4), fliplr(log10(idx4))], ...
      [log10(max(m4(idx4)+s4(idx4), 1e-10)), fliplr(log10(max(m4(idx4)-s4(idx4), 1e-10)))], ...
      'cyan', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(log10(idx4), log10(m4(idx4)), ...
    'LineStyle','-', 'Color','cyan', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-5}~\mathrm{Adaptive}~\eta_i=\frac{1}{2+25000\mathrm{Tr}(P_{\hat{K}_i})}$');

% --- Sigma_w = 1e-3, adaptive step size ---
patch([log10(idx5), fliplr(log10(idx5))], ...
      [log10(max(m5(idx5)+s5(idx5), 1e-10)), fliplr(log10(max(m5(idx5)-s5(idx5), 1e-10)))], ...
      'magenta', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(log10(idx5), log10(m5(idx5)), ...
    'LineStyle','--', 'Color','magenta', 'LineWidth',1.5, ...
    'DisplayName', '$\Sigma_w=10^{-3}~\mathrm{Adaptive}~\eta_i=\frac{1}{2+25000\mathrm{Tr}(P_{\hat{K}_i})}$');

% --- Sigma_w = 1e-1, adaptive step size (currently disabled) ---
% idx6 = 1:step:numel(m6);
% patch([log10(idx6), fliplr(log10(idx6))], ...
%       [log10(max(m6(idx6)+s6(idx6), 1e-10)), fliplr(log10(max(m6(idx6)-s6(idx6), 1e-10)))], ...
%       'magenta', 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
% plot(log10(idx6), log10(m6(idx6)), ...
%     'LineStyle','-', 'Color','magenta', 'LineWidth',1.5, ...
%     'DisplayName', '$\Sigma_w=10^{-1}~Adaptive~\eta_i=\frac{1}{2+25000\mathrm{Tr}(P_{\hat{K}_i})}$');

% Manually set axis ticks to display original (non-log) values, since the
% data itself is plotted on a log10 scale
ax = gca;
ax.XAxis.TickValues = log10([1, 10, 100, 1000, 10000, 100000, 500000]);
ax.XAxis.TickLabels = {'1','10','100','10^3','10^4','10^5','5\times10^5'};
ax.YAxis.TickValues = log10([0.001, 0.01, 0.1, 1, 10, 100]);
ax.YAxis.TickLabels = {'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'};
ylim([log10(0.05), log10(1.1)]);
xlim([log10(10),log10(300000)])

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