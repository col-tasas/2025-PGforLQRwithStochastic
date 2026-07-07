%% ========================================================================
%  不同 noise level 下比较 baseline vs. non-baseline 梯度估计准确性
%  公平对比: 总样本预算相同 (n=500)
%    - No-baseline:      n = 500 条扰动轨迹
%    - With baseline:    n = 400 条扰动轨迹 + 100 条 baseline 轨迹 = 500 总数
%  Ground truth: LQR 解析策略梯度
%  ========================================================================
clear; clc; close all;

%% ---------------- 系统设置 ----------------
A=[1,-1.13,-0.65,-0.807,1.59; 0,0.77,0.32,-0.98,-2.97;0,0.12,0.02,0.0,-0.36;0,0.01,0.01,-0.03,-0.04;0,0.14,-0.09,0.29,0.76];
B=[89.20,-50.17,1.13,-19.35;5.22,6.36,0.23,-0.32;-9.47,5.93,-0.12,0.99;-0.32,0.32,-0.01,-0.01;-4.53,3.21,-0.14,0.09];
nx=size(A,2);
nu=size(B,2);
Q=eye(5);
R=eye(4);

[~,K_initial]=idare(A,B,40*Q,R);
K_initial=-K_initial;

disp('闭环特征值 (应全部在单位圆内):');
disp(eig(A+B*K_initial));

sqrtSigma_0 = sqrt(1e-6);
Sigma_0     = 1e-6*eye(nx);

l = 20;     % rollout 长度 (固定)
r = 0.01;   % 扰动半径 (固定)

n_perturb_baseline = 400;   % 带 baseline 方法: 扰动轨迹数
n_baseline_samples = 100;   % 带 baseline 方法: baseline 估计轨迹数 (400+100=500)
n_nobaseline       = 500;   % 不带 baseline 方法: 扰动轨迹数 (总预算同为500)

%% ---------------- 噪声等级扫描 ----------------
noise_levels = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];  % Sigma_w = level * I
nTrials = 100;

err_baseline = zeros(length(noise_levels), nTrials);
err_nobase   = zeros(length(noise_levels), nTrials);

fprintf('=== 误差 vs. noise level (总样本预算 n=500 固定) ===\n');
for i = 1:length(noise_levels)
    level = noise_levels(i);
    Sigma_w     = level*eye(nx);
    sqrtSigma_w = sqrt(level);

    % ground truth 梯度 (随噪声等级变化, 因为 Sigma_K 依赖 Sigma_w)
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

%% ---------------- 画图: Boxplot 展示方差 (y轴 log scale) ----------------
% 直接用原始相对误差画 boxplot, 然后把 y 轴设成 log scale (而不是手动
% log10 变换数据) —— y 轴刻度标签直接显示真实的相对误差数值, 是真正的
% "loglog" 风格 (x轴noise level本身也可以设成log刻度).
%
% 区分两种方法的方式:
%   - With baseline:    箱子填色 (蓝色), outlier 用 'x' 标记
%   - Without baseline: 箱子不填色 (只有黑色边框), outlier 用 'o' 标记
nLevels = length(noise_levels);

data  = [];
group = [];

groupWidth = 0.01;   % 同一 noise level 内, 左右两个箱子的间距
for i = 1:nLevels
    data  = [data, err_baseline(i,:), err_nobase(i,:)];
    group = [group, ...
        repmat(2*i-1, 1, nTrials), ...   % baseline 用奇数编号 (左)
        repmat(2*i,   1, nTrials)];      % no-baseline 用偶数编号 (右)
end

figure('Position',[100 100 900 520]);
positions = zeros(1, 2*nLevels);
for i = 1:nLevels
    positions(2*i-1) = i - groupWidth/2;   % 左
    positions(2*i)   = i + groupWidth/2;   % 右
end

boxplot(data, group, 'positions', positions, 'widths', 0.28, ...
    'Colors', 'k', 'Symbol','');

set(gca, 'YScale', 'log');   % 关键: y 轴 log scale, 数据本身不做log变换

% --- 上色 (只给 with-baseline 填色) ---
% findobj 返回顺序是"倒序" (最后一组画在最前面),
% boxObj(j) 对应真实 group 编号是 (n - j + 1).
boxObj = findobj(gca,'Tag','Box');
nBoxes = length(boxObj);

fillColor = [0.00 0.45 0.74];  % with-baseline 填色蓝

for j = 1:nBoxes
    groupNum = nBoxes - j + 1;
    isBaseline = mod(groupNum,2) == 1;   % 奇数(左) = with baseline
    if isBaseline
        patch(get(boxObj(j),'XData'), get(boxObj(j),'YData'), fillColor, 'FaceAlpha',0.6);
        hold on
    end
    % without baseline: 不填色, 保留boxplot默认画的黑色边框即可
end

% --- 连接每种方法各个 noise level 箱子的中位数 ---
medianObj = findobj(gca,'Tag','Median');
nMedians  = length(medianObj);
medianX_baseline = zeros(1,nLevels);
medianY_baseline = zeros(1,nLevels);
medianX_nobase   = zeros(1,nLevels);
medianY_nobase   = zeros(1,nLevels);

for j = 1:nMedians
    groupNum = nMedians - j + 1;         % 该 median 对应的真实 group 编号
    xd = get(medianObj(j),'XData');
    yd = get(medianObj(j),'YData');
    xMid = mean(xd);
    yMid = yd(1);   % median 线是水平的, 两端 y 值相同
    levelIdx = ceil(groupNum/2);          % group 1,2 -> level1; group 3,4 -> level2; ...
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
%ylabel('$\frac{\lVert \hat{\nabla} C(\hat{K}_0) - \nabla C(\hat{K}_0)\rVert_F] {\lVert\nabla C(\hat{K}_0)\rVert_F}$','interpreter','latex','FontSize',12);
%ylabel('$\frac{\lVert \hat{\nabla} C(\hat{K}_0) - \nabla C(\hat{K}_0)\rVert_F}{\lVert\nabla C(\hat{K}_0)\rVert_F}$','interpreter','latex','FontSize',12)
ylabel('$\frac{\| \hat{\nabla} C(\hat{K}_0) - \nabla C(\hat{K}_0)\|_F}{\|\nabla C(\hat{K}_0)\|_F}$', 'Interpreter','latex','FontSize',12)
%title(sprintf('梯度估计误差分布 vs. Noise Level  (r=%.3f, l=%d, 总预算=500, %d trials)', r, l, nTrials));
grid on; hold on;

% 图例
h1 = patch(NaN,NaN,fillColor,'FaceAlpha',0.6);
h2 = patch(NaN,NaN,'w','EdgeColor','k');
legend([h2 h1], { sprintf('Algorithm 1 : $n=500$')
    sprintf('Algorithm 3 (VR): $n_b=400,n_v=100$')},'Interpreter','latex','FontSize',12, 'Location','best');

%% ========================================================================
%  函数定义
%  ========================================================================

function nablaK_true = trueGradient(A,B,Q,R,K,Sigma_0,Sigma_w)
% 解析 LQR 策略梯度
    AK = A + B*K;
    P     = dlyap(AK', Q + K'*R*K);
    Sigma = dlyap(AK, Sigma_0 + Sigma_w);
    nablaK_true = 2*((R + B'*P*B)*K + B'*P*A) * Sigma;
end

function nablaK = estimationBS_vec(kcurrent,l,r,n,n_baseline,A,B,Q,R,sqrtSigma_w,sqrtSigma_0,nx,nu)
% 带 baseline 的零阶梯度估计 (方差缩减)
% n:          扰动轨迹数
% n_baseline: 用于估计 baseline 的轨迹数
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
% 不带 baseline 的零阶梯度估计
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
% 固定初始状态 x_0 (通常为0), 向量化 n 条轨迹算 baseline cost
    x  = repmat(x_0, [1, n]);
    ck = zeros(1, n);
    for p = 1:l
        u  = kcurrent * x;
        ck = ck + sum(x.*(Q*x),1) + sum(u.*(R*u),1);
        x  = A*x + B*u + sqrtSigma_w*randn(nx,n);
    end
    b = mean(ck) / l;
end
