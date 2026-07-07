%% Define the true system and the weighting matrices
A=[1,-1.13,-0.65,-0.807,1.59; 0,0.77,0.32,-0.98,-2.97;0,0.12,0.02,0.0,-0.36;0,0.01,0.01,-0.03,-0.04;0,0.14,-0.09,0.29,0.76]; % system matrix A
B=[89.20,-50.17,1.13,-19.35;5.22,6.36,0.23,-0.32;-9.47,5.93,-0.12,0.99;-0.32,0.32,-0.01,-0.01;-4.53,3.21,-0.14,0.09]; % input matrix B
nx=size(A,2); % number of states
nu=size(B,2); % number of control inputs
Q=eye(5); % weighting matrix Q
R=eye(4); % weighting matrix R

% Initial stabilizing gain (from DARE with inflated state cost)
[Pstar1,K_initial]=idare(A,B,40*Q,R);
% Optimal LQR cost-to-go (ground truth for comparison)
[Pstar,~]=idare(A,B,Q,R);
K_initial=-K_initial;
eig(A+B*K_initial)   % check closed-loop stability (should be inside unit circle)

Sigma_w=10^(-5)*eye(nx);       % process noise covariance
Sigma_0=10^(-6)*eye(nx);       % initial state covariance
% Pre-compute the standard deviations (Sigma_w and Sigma_0 are scalar
% multiples of eye(nx), so sampling via scalar*randn avoids repeated
% calls to mvnrnd with a full covariance matrix)
sqrtSigma_w = sqrt(10^(-5));
sqrtSigma_0 = sqrt(10^(-6));

iteration=500000;   % number of PG update steps
numup=100;              % number of independent trials
ck=zeros(numup,iteration);  % normalized sub-optimality gap (only row 1 is used here, since numup=1)

n=500;    % number of perturbed rollouts used for gradient/covariance estimation
l=20;     % rollout length
r=0.01;   % perturbation radius

pk=dlyap((A+B*K_initial)',Q+K_initial'*R*K_initial);
eta=a / (b + c*trace(pk));
%% Natural Policy Gradient (NPG)
for num=1:numup
    kcurrentBS=K_initial;
    for i=1:iteration
        [nablaK,SigmaK]=estimationBS_vec(kcurrentBS,l,r,n,A,B,Q,R,sqrtSigma_w,sqrtSigma_0,nx,nu);
        pk=dlyap((A+B*kcurrentBS)',Q+kcurrentBS'*R*kcurrentBS);
        ck(num,i)=(trace(pk)-trace(Pstar))/trace(Pstar);
        disp(['i = ', num2str(i), ', ck = ', num2str(ck(num,i))]);
        % ---- Natural Policy Gradient update ----
        kcurrentBS = kcurrentBS - eta*nablaK/SigmaK;
    end
end
save("NPG3NonAda.mat")

%% Required functions
% Vectorized gradient/covariance estimation: all n samples are simulated
% simultaneously; only the time loop over p = 1:l remains.
function [nablaK,SigmaK]=estimationBS_vec(kcurrent,l,r,n,A,B,Q,R,sqrtSigma_w,sqrtSigma_0,nx,nu)
    % true kernel P_K (diagnostic only, not used in the gradient update)
    P=dlyap((A+B*kcurrent)',Q+kcurrent'*R*kcurrent);
    % true covariance Sigma_K (diagnostic only)
    % Pw=dlyap((A+B*kcurrent),Sigma_w);
    % C=2*((R+B'*P*B)*kcurrent+B'*P*A)*Pw;

    % ---- Vectorized sampling: generate all n samples at once ----
    x = sqrtSigma_0 * randn(nx, n);              % nx x n, initial state (each column is one sample)
    Xpert = randn(nu*nx, n);
    Xpert = Xpert ./ vecnorm(Xpert, 2, 1);        % normalize each column to unit norm (uniform sampling on a sphere)
    Xpert = r * Xpert;                             % scale to radius r
    Xpert3 = reshape(Xpert, [nu, nx, n]);          % nu x nx x n
    Kpert3 = repmat(kcurrent, [1,1,n]) + Xpert3;   % nu x nx x n, each page is one perturbed gain K_j
    ck = zeros(1, n);             % accumulated cost for each sample
    xk = zeros(nx, nx, n);        % accumulated x*x' for each sample (used for the SigmaK estimate)
    for p = 1:l
        x3 = reshape(x, [nx, 1, n]);                  % nx x 1 x n
        u3 = pagemtimes(Kpert3, x3);                   % nu x 1 x n, batched matrix multiply K_j*x_j
        u = reshape(u3, [nu, n]);                       % nu x n
        Qx = Q * x;        % nx x n
        Ru = R * u;        % nu x n
        ck = ck + sum(x .* Qx, 1) + sum(u .* Ru, 1);    % 1 x n, equivalent to x'(Q+K'RK)x
        xk = xk + reshape(x, [nx,1,n]) .* reshape(x, [1,nx,n]);  % nx x nx x n
        noise = sqrtSigma_w * randn(nx, n);
        x = A*x + B*u + noise;                          % nx x n, state update (A,B shared, broadcast across columns)
    end
    % Gradient estimate: (nx*nu/r^2) * (1/l) * sum_j ck(j) * X_j, then
    % averaged over n
    scale = (nx*nu/(r^2)) * (1/l) * ck;                 % 1 x n
    weighted = reshape(scale, [1,1,n]) .* Xpert3;        % nu x nx x n
    nablaK = sum(weighted, 3) / n;                       % nu x nx
    % errorgradient = norm(nablaK - C);
    SigmaK = sum(xk, 3) / (n*l);
    % errorcovariance = norm(SigmaK - Pw);
end