%% ========================================================================
%  Compare baseline vs. non-baseline zeroth-order gradient estimation
%  accuracy across different noise levels.
%  Fair comparison: identical total sample budget (n = 500)
%    - No-baseline:   n = 500 perturbed rollouts
%    - With baseline: n = 400 perturbed rollouts + 100 baseline rollouts
%                      = 500 total
%  Ground truth: analytical LQR policy gradient
%  ========================================================================
clear; clc; close all;

%% ---------------- System setup ----------------
A=[1,-1.13,-0.65,-0.807,1.59; 0,0.77,0.32,-0.98,-2.97;0,0.12,0.02,0.0,-0.36;0,0.01,0.01,-0.03,-0.04;0,0.14,-0.09,0.29,0.76];
B=[89.20,-50.17,1.13,-19.35;5.22,6.36,0.23,-0.32;-9.47,5.93,-0.12,0.99;-0.32,0.32,-0.01,-0.01;-4.53,3.21,-0.14,0.09];
nx=size(A,2);
nu=size(B,2);
Q=eye(5);
R=eye(4);

% Initial stabilizing gain (from DARE with inflated state cost)
[~,K_initial]=idare(A,B,40*Q,R);
K_initial=-K_initial;

disp('Closed-loop eigenvalues (should all lie inside the unit circle):');
disp(eig(A+B*K_initial));

sqrtSigma_0 = sqrt(1e-6);
Sigma_0     = 1e-6*eye(nx);

l = 20;     % rollout length (fixed)
r = 0.01;   % perturbation radius (fixed)

n_perturb_baseline = 400;   % with-baseline method: number of perturbed rollouts
n_baseline_samples = 100;   % with-baseline method: number of baseline-estimation rollouts (400+100=500)
n_nobaseline       = 500;   % no-baseline method: number of perturbed rollouts (same total budget of 500)

%% ---------------- Sweep over noise levels ----------------
noise_levels = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];  % Sigma_w = level * I
nTrials = 100;

err_baseline = zeros(length(noise_levels), nTrials);
err_nobase   = zeros(length(noise_levels), nTrials);

fprintf('=== Error vs. noise level (fixed total sample budget n=500) ===\n');
for i = 1:length(noise_levels)
    level = noise_levels(i);
    Sigma_w     = level*eye(nx);
    sqrtSigma_w = sqrt(level);

    % ground-truth gradient (depends on the noise level, since Sigma_K
    % depends on Sigma_w)
    nablaK_true = trueGradient(A,B,Q,R,K_initial,Sigma_0,Sigma_w);
    trueNorm = norm(nablaK_true,'fro');

    for t = 1:nTrials
        g1 = estimationBS_vec(K_initial,l,r,n_perturb_baseline,n_baseline_samples, ...
                A,B,Q,R,sqrtSigma_w,sqrtSigma_0,nx,nu);
        g2 = estimationBS_vec_nobaseline(K_initial,l,r,n_nobaseline, ...
                A,B,Q,R,sqrtSigma_w,sqrtSigma_0,nx,nu);

        err_baseline(i,t) = norm(g1 - nablaK_true) / trueNorm;
        err_nobase(i,t)   = norm(g2 - nablaK_true) / trueNorm;
    end
    fprintf('noise=%.1e  |  baseline err = %.4e  |  no-baseline err = %.4e\n', ...
        level, mean(err_baseline(i,:)), mean(err_nobase(i,:)));
end

mean_err_baseline = mean(err_baseline,2);
mean_err_nobase   = mean(err_nobase,2);
std_err_baseline  = std(err_baseline,0,2)/sqrt(nTrials);
std_err_nobase    = std(err_nobase,0,2)/sqrt(nTrials);

%% ---------------- Plot: boxplot showing variance (log-scale y-axis) ----------------
% Boxplot is drawn directly on the raw relative errors, then the y-axis
% is set to log scale (rather than manually log10-transforming the data)
% so the axis tick labels show the true relative-error values -- this
% gives a genuine "loglog" style plot (the x-axis noise level can also be
% set to a log scale).
%
% The two methods are distinguished visually as follows:
%   - With baseline:    box is filled (blue), outliers marked with 'x'
%   - Without baseline: box is unfilled (black outline only), outliers marked with 'o'
nLevels = length(noise_levels);

data  = [];
group = [];

groupWidth = 0.01;   % spacing between the two boxes within the same noise level
for i = 1:nLevels
    data  = [data, err_baseline(i,:), err_nobase(i,:)];
    group = [group, ...
        repmat(2*i-1, 1, nTrials), ...   % baseline uses odd indices (left box)
        repmat(2*i,   1, nTrials)];      % no-baseline uses even indices (right box)
end

figure('Position',[100 100 900 520]);
positions = zeros(1, 2*nLevels);
for i = 1:nLevels
    positions(2*i-1) = i - groupWidth/2;   % left
    positions(2*i)   = i + groupWidth/2;   % right
end

boxplot(data, group, 'positions', positions, 'widths', 0.28, ...
    'Colors', 'k', 'Symbol','');

set(gca, 'YScale', 'log');   % key: log-scale y-axis, data itself is not log-transformed

% --- Color fill (only the with-baseline boxes get filled) ---
% findobj returns objects in reverse order (the last group drawn appears
% first), so boxObj(j) corresponds to the actual group number (n - j + 1).
boxObj = findobj(gca,'Tag','Box');
nBoxes = length(boxObj);

fillColor = [0.00 0.45 0.74];  % fill color (blue) for with-baseline boxes

for j = 1:nBoxes
    groupNum = nBoxes - j + 1;
    isBaseline = mod(groupNum,2) == 1;   % odd (left) = with baseline
    if isBaseline
        patch(get(boxObj(j),'XData'), get(boxObj(j),'YData'), fillColor, 'FaceAlpha',0.6);
        hold on
    end
    % without baseline: leave unfilled, keep boxplot's default black outline
end

% --- Connect the median of each noise level's box for each method ---
medianObj = findobj(gca,'Tag','Median');
nMedians  = length(medianObj);
medianX_baseline = zeros(1,nLevels);
medianY_baseline = zeros(1,nLevels);
medianX_nobase   = zeros(1,nLevels);
medianY_nobase   = zeros(1,nLevels);

for j = 1:nMedians
    groupNum = nMedians - j + 1;         % actual group number for this median
    xd = get(medianObj(j),'XData');
    yd = get(medianObj(j),'YData');
    xMid = mean(xd);
    yMid = yd(1);   % the median line is horizontal, so both endpoints share the same y value
    levelIdx = ceil(groupNum/2);          % groups 1,2 -> level 1; groups 3,4 -> level 2; ...
    if mod(groupNum,2) == 1
        medianX_baseline(levelIdx) = xMid;
        medianY_baseline(levelIdx) = yMid;
    else
        medianX_nobase(levelIdx) = xMid;
        medianY_nobase(levelIdx) = yMid;
    end
end

plot(medianX_baseline, medianY_baseline, '-', 'Color', fillColor, 'LineWidth', 1.5);
hold on
plot(medianX_nobase,   medianY_nobase,   '--', 'Color', 'k',       'LineWidth', 1.5);
hold on
set(gca, 'XTick', 1:nLevels, 'XTickLabel', arrayfun(@(x) sprintf('%.0e',x), noise_levels, 'UniformOutput', false));
xlabel('Noise level  $\Sigma_w = \#\times I_5$','interpreter','latex','FontSize',12);
ylabel('$\frac{\| \hat{\nabla} C(\hat{K}_0) - \nabla C(\hat{K}_0)\|_F}{\|\nabla C(\hat{K}_0)\|_F}$', 'Interpreter','latex','FontSize',12)
grid on; hold on;

% Legend
h1 = patch(NaN,NaN,fillColor,'FaceAlpha',0.6);
h2 = patch(NaN,NaN,'w','EdgeColor','k');
legend([h2 h1], { sprintf('Algorithm 1 : $n=500$')
    sprintf('Algorithm 3 (VR): $n_b=400,n_v=100$')},'Interpreter','latex','FontSize',12, 'Location','best');

%% ========================================================================
%  Function definitions
%  ========================================================================

function nablaK_true = trueGradient(A,B,Q,R,K,Sigma_0,Sigma_w)
% Analytical LQR policy gradient
    AK = A + B*K;
    P     = dlyap(AK', Q + K'*R*K);
    Sigma = dlyap(AK, Sigma_0 + Sigma_w);
    nablaK_true = 2*((R + B'*P*B)*K + B'*P*A) * Sigma;
end

function nablaK = estimationBS_vec(kcurrent,l,r,n,n_baseline,A,B,Q,R,sqrtSigma_w,sqrtSigma_0,nx,nu)
% Zeroth-order gradient estimate with baseline (variance reduction)
% n:          number of perturbed rollouts
% n_baseline: number of rollouts used to estimate the baseline
    bfunction = estimationBaseline(zeros(nx,1),kcurrent,l,n_baseline,A,B,Q,R,sqrtSigma_w,nx,nu);

    x = sqrtSigma_0 * randn(nx, n);
    Xpert = randn(nu*nx, n);
    Xpert = Xpert ./ vecnorm(Xpert, 2, 1);
    Xpert = r * Xpert;
    Xpert3 = reshape(Xpert, [nu, nx, n]);
    Kpert3 = repmat(kcurrent, [1,1,n]) + Xpert3;
    ck = zeros(1, n);
    for p = 1:l
        x3 = reshape(x, [nx, 1, n]);
        u  = reshape(pagemtimes(Kpert3, x3), [nu, n]);
        ck = ck + sum(x.*(Q*x),1) + sum(u.*(R*u),1);
        x  = A*x + B*u + sqrtSigma_w*randn(nx,n);
    end
    ck = ck / l;
    scale    = (nx*nu/(r^2)) * (ck - bfunction);
    weighted = reshape(scale, [1,1,n]) .* Xpert3;
    nablaK   = sum(weighted, 3) / n;
end

function nablaK = estimationBS_vec_nobaseline(kcurrent,l,r,n,A,B,Q,R,sqrtSigma_w,sqrtSigma_0,nx,nu)
% Zeroth-order gradient estimate without baseline
    x = sqrtSigma_0 * randn(nx, n);
    Xpert = randn(nu*nx, n);
    Xpert = Xpert ./ vecnorm(Xpert, 2, 1);
    Xpert = r * Xpert;
    Xpert3 = reshape(Xpert, [nu, nx, n]);
    Kpert3 = repmat(kcurrent, [1,1,n]) + Xpert3;
    ck = zeros(1, n);
    for p = 1:l
        x3 = reshape(x, [nx, 1, n]);
        u3 = pagemtimes(Kpert3, x3);
        u  = reshape(u3, [nu, n]);
        ck = ck + sum(x .* (Q*x), 1) + sum(u .* (R*u), 1);
        x  = A*x + B*u + sqrtSigma_w*randn(nx, n);
    end
    scale    = (nx*nu/(r^2)) * (1/l) * ck;
    weighted = reshape(scale, [1,1,n]) .* Xpert3;
    nablaK   = sum(weighted, 3) / n;
end

function b = estimationBaseline(x_0,kcurrent,l,n,A,B,Q,R,sqrtSigma_w,nx,nu)
% Fixed initial state x_0 (typically zero); vectorized over n rollouts to
% compute the baseline cost
    x  = repmat(x_0, [1, n]);
    ck = zeros(1, n);
    for p = 1:l
        u  = kcurrent * x;
        ck = ck + sum(x.*(Q*x),1) + sum(u.*(R*u),1);
        x  = A*x + B*u + sqrtSigma_w*randn(nx,n);
    end
    b = mean(ck) / l;
end
