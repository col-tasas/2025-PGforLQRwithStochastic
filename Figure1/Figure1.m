%% Policy Gradient for LQR: Effect of Gradient Noise and Step Size
% This script studies zeroth/first-order policy gradient (PG) updates for
% the discrete-time LQR problem under different levels of injected
% gradient noise and different step sizes. It runs four experiments:
%   1. Noise-free PG (single trajectory, deterministic)
%   2. PG with low noise variance, large step size
%   3. PG with high noise variance, large step size
%   4. PG with high noise variance, small step size
% Each experiment saves its results to a .mat file for later analysis.

%% Problem setup
A = [1, -1.13, -0.65, -0.807, 1.59;
     0, 0.77, 0.32, -0.98, -2.97;
     0, 0.12, 0.02, 0.0, -0.36;
     0, 0.01, 0.01, -0.03, -0.04;
     0, 0.14, -0.09, 0.29, 0.76];

B = [89.20, -50.17, 1.13, -19.35;
     5.22, 6.36, 0.23, -0.32;
     -9.47, 5.93, -0.12, 0.99;
     -0.32, 0.32, -0.01, -0.01;
     -4.53, 3.21, -0.14, 0.09];

nx = size(A, 2);   % state dimension
nu = size(B, 2);   % input dimension

Q = eye(5);         % state cost weight
R = eye(4);         % input cost weight
Sigma_w = 10^(-1) * eye(nx);  % process noise covariance used in Lyapunov weighting

% Initial stabilizing gain: solve DARE with an inflated state cost (40*Q)
% to get a more conservative (better-damped) starting policy
[Pstar1, K_initial] = idare(A, B, 40*Q, R);
K_initial = -K_initial;

% Optimal LQR gain and cost-to-go (ground truth for comparison)
[Pstar, Kstar] = idare(A, B, Q, R);
Kstar = -Kstar;
trPstar = trace(Pstar);

%% Experiment parameters
iteration = 500000;      % number of PG update steps
number = 100;             % number of parallel sample trajectories (Monte Carlo runs)
step_high = 0.5 * 10^(-3);  % large step size
step_low = 0.2 * 10^(-3);   % small step size
noisevar_high = 5;        % high gradient noise variance
noisevar_low = 1;         % low gradient noise variance

%% Experiment 1: Noise-free PG (single trajectory, no randomness, no averaging needed)
disp('=== Running noise-free PG ===');
CKPGV = zeros(1, iteration);
K = K_initial;
for i = 1:iteration
    P  = dlyap((A + B*K)', Q + K'*R*K);   % cost-to-go Lyapunov solution
    Pw = dlyap((A + B*K), Sigma_w);        % state covariance Lyapunov solution
    CKPGV(i) = trace(P - Pstar) / trPstar; % normalized sub-optimality gap
    K = K - step_high * ((R + B'*P*B)*K + B'*P*A) * Pw;  % exact PG update
    if mod(i, 1000) == 0
        disp(['i=', num2str(i), ', ck=', num2str(CKPGV(i))]);
    end
end
save('dataPG_noisefree.mat', 'CKPGV', 'step_high');
disp('Noise-free done and saved.');

%% Experiment 2: low noise variance, large step size
disp('=== Running PG2: noisevar_low, step_high ===');
ck_mean = zeros(1, iteration);
ck_std  = zeros(1, iteration);
KPG = repmat(K_initial, [1, 1, number]);
for i = 1:iteration
    KPG_i = KPG;
    ck = zeros(number, 1);
    parfor j = 1:number
        K1  = KPG_i(:, :, j);
        P1  = dlyap((A + B*K1)', Q + K1'*R*K1);
        Pw1 = dlyap((A + B*K1), Sigma_w);
        ck(j) = trace(P1 - Pstar) / trPstar;
        % noisy PG update: add zero-mean Gaussian noise to the gradient
        KPG_i(:, :, j) = K1 - step_high * (((R + B'*P1*B)*K1 + B'*P1*A) * Pw1 ...
            + normrnd(0, noisevar_low, [nu, nx]));
    end
    KPG = KPG_i;
    ck_mean(i) = mean(ck);
    ck_std(i)  = std(ck);
    if mod(i, 1000) == 0
        disp(['i=', num2str(i), ', mean=', num2str(ck_mean(i))]);
        save('dataPG_noiselow_stephigh.mat', 'ck_mean', 'ck_std', 'step_high', 'noisevar_low');
    end
end
save('dataPG_noiselow_stephigh.mat', 'ck_mean', 'ck_std', 'step_high', 'noisevar_low');
disp('PG2 done and saved.');

%% Experiment 3: high noise variance, large step size
disp('=== Running PG3: noisevar_high, step_high ===');
ck_mean = zeros(1, iteration);
ck_std  = zeros(1, iteration);
KPG = repmat(K_initial, [1, 1, number]);
for i = 1:iteration
    KPG_i = KPG;
    ck = zeros(number, 1);
    parfor j = 1:number
        K1  = KPG_i(:, :, j);
        P1  = dlyap((A + B*K1)', Q + K1'*R*K1);
        Pw1 = dlyap((A + B*K1), Sigma_w);
        ck(j) = trace(P1 - Pstar) / trPstar;
        if ck(j) < 0
            % guard against numerical instability producing a negative
            % (i.e. invalid) sub-optimality gap; flag with a large value
            ck(j) = 100;
        end
        KPG_i(:, :, j) = K1 - step_high * (((R + B'*P1*B)*K1 + B'*P1*A) * Pw1 ...
            + normrnd(0, noisevar_high, [nu, nx]));
    end
    KPG = KPG_i;
    ck_mean(i) = mean(ck);
    ck_std(i)  = std(ck);
    if mod(i, 1000) == 0
        disp(['i=', num2str(i), ', mean=', num2str(ck_mean(i))]);
        save('dataPG_noisehigh_stephigh.mat', 'ck_mean', 'ck_std', 'step_high', 'noisevar_high');
    end
end
save('dataPG_noisehigh_stephigh.mat', 'ck_mean', 'ck_std', 'step_high', 'noisevar_high');
disp('PG3 done and saved.');

%% Experiment 4: high noise variance, small step size
disp('=== Running PG4: noisevar_high, step_low ===');
ck_mean = zeros(1, iteration);
ck_std  = zeros(1, iteration);
KPG = repmat(K_initial, [1, 1, number]);
for i = 1:iteration
    KPG_i = KPG;
    ck = zeros(number, 1);
    parfor j = 1:number
        K1  = KPG_i(:, :, j);
        P1  = dlyap((A + B*K1)', Q + K1'*R*K1);
        Pw1 = dlyap((A + B*K1), Sigma_w);
        ck(j) = trace(P1 - Pstar) / trPstar;
        KPG_i(:, :, j) = K1 - step_low * (((R + B'*P1*B)*K1 + B'*P1*A) * Pw1 ...
            + normrnd(0, noisevar_high, [nu, nx]));
    end
    KPG = KPG_i;
    ck_mean(i) = mean(ck);
    ck_std(i)  = std(ck);
    if mod(i, 1000) == 0
        disp(['i=', num2str(i), ', mean=', num2str(ck_mean(i))]);
        save('dataPG_noisehigh_steplow.mat', 'ck_mean', 'ck_std', 'step_low', 'noisevar_high');
    end
end
save('dataPG_noisehigh_steplow.mat', 'ck_mean', 'ck_std', 'step_low', 'noisevar_high');
disp('PG4 done and saved.');

disp('All done.');
