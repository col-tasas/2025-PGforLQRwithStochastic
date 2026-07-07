%% Define the true system and the weighting matrices
A=[1,-1.13,-0.65,-0.807,1.59; 0,0.77,0.32,-0.98,-2.97;0,0.12,0.02,0.0,-0.36;0,0.01,0.01,-0.03,-0.04;0,0.14,-0.09,0.29,0.76];
B=[89.20,-50.17,1.13,-19.35;5.22,6.36,0.23,-0.32;-9.47,5.93,-0.12,0.99;-0.32,0.32,-0.01,-0.01;-4.53,3.21,-0.14,0.09];
nx=size(A,2);
nu=size(B,2);
Q=eye(5);
R=eye(4);

% Initial stabilizing gain (from DARE with inflated state cost)
[Pstar1,K_initial]=idare(A,B,40*Q,R);
% Optimal LQR cost-to-go (ground truth for comparison)
[Pstar,~]=idare(A,B,Q,R);
K_initial=-K_initial;
eig(A+B*K_initial)   % check closed-loop stability (should be inside unit circle)

sqrtSigma_w = sqrt(10^(-3));   % process noise std
sqrtSigma_0 = sqrt(10^(-6));   % initial state std

iteration = 3500000;   % number of PG update steps per trial
numup = 100;              % number of independent trials
ck = zeros(numup, iteration);  % normalized sub-optimality gap, per trial/iteration

n   = 400;      % number of perturbed rollouts used for gradient estimation
l   = 20;       % rollout length
r   = 0.01;     % perturbation radius
eta = 0.0008;   % step size

%% Zeroth-order policy gradient with baseline (variance reduction)
for num = 1:numup
    kcurrentBS = K_initial;
    for i = 1:iteration
        try
            nablaK = estimationBS_vec(kcurrentBS,l,r,n,A,B,Q,R,sqrtSigma_w,sqrtSigma_0,nx,nu);
            pk = dlyap((A+B*kcurrentBS)', Q+kcurrentBS'*R*kcurrentBS);
            ck(num,i) = (trace(pk)-trace(Pstar))/trace(Pstar);
            if mod(i,1000)==0
                disp(['Trial ', num2str(num), '/', num2str(numup), ...
                    ', i = ', num2str(i), ...
                    ', ck = ', num2str(ck(num,i))]);
            end
            kcurrentBS = kcurrentBS - eta*nablaK;
        catch
            % stop this trial early if the update becomes numerically
            % invalid (e.g. Lyapunov equation fails to solve because the
            % gain destabilized the closed loop)
            break
        end
    end
    save('PG_VR3.mat', 'ck');
    disp(['Trial ', num2str(num), ' saved.']);
end
save('PG_VR3.mat', 'ck');

%% Required functions
function nablaK = estimationBS_vec(kcurrent,l,r,n,A,B,Q,R,sqrtSigma_w,sqrtSigma_0,nx,nu)
% Zeroth-order gradient estimate with baseline (variance reduction),
% vectorized over n perturbed rollouts

    % baseline: fixed x_0 = 0, rolled out with the unperturbed gain K
    bfunction = estimationBaseline(zeros(nx,1),kcurrent,l,100,A,B,Q,R,sqrtSigma_w,nx,nu);

    % initial state, nx x n
    x = sqrtSigma_0 * randn(nx, n);

    % perturbation directions: uniform sampling on a sphere, scaled to
    % radius r
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

function b = estimationBaseline(x_0,kcurrent,l,n,A,B,Q,R,sqrtSigma_w,nx,nu)
% Fixed initial state x_0, vectorized over n rollouts
    x  = repmat(x_0, [1, n]);   % nx x n
    ck = zeros(1, n);
    for p = 1:l
        u  = kcurrent * x;      % nu x n
        ck = ck + sum(x.*(Q*x),1) + sum(u.*(R*u),1);
        x  = A*x + B*u + sqrtSigma_w*randn(nx,n);
    end
    b = mean(ck) / l;
end