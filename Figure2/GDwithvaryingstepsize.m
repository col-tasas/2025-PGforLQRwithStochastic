%% Parameters 
% % define the true system and the weighting matrix 
eta_BS_4=[];
eta_BS_2=[1];
eta_BS_0=[];
eta_BS=[eta_BS_4,eta_BS_2,eta_BS_0];
episode=100;
varying=zeros(1,episode);
for i=1:15
    varying(1,i)=1.08^i;
end
% for i=16:episode
%     varying(1,i)=varying(1,10);
% end 
AmpNoise_forloop_BS=[0.0001*ones(1,numel(eta_BS_4)),0.01*ones(1,numel(eta_BS_2)),1*ones(1,numel(eta_BS_0))];
Sigma_0=0.0001*eye(3);
%%
A=[1.01,0.01,0;0.01,1.01,0.01;0,0.01,1.01];% system matrix A
B=1*eye(3); % input matrix B
nx=size(A,2); % number of states
nu=size(B,2); % number of control inputs
Q=0.001*eye(3); % weighting matrix Q
R=1*eye(3); % weighting matrix R
[Pstar1,K_initial]=idare(A,B,50*Q,R);
K_initial=-K_initial;

num_BS=numel(AmpNoise_forloop_BS);
%% Algorithm
% K optimal
[Pstar,Kstar]=idare(A,B,Q,R);
Kstar=-Kstar;
errorgradientBS=zeros(episode,num_BS);
errorcovarianceBS=zeros(episode,num_BS);
khistoryBS=zeros(nx*nu,episode,num_BS);
SigmahistoryBS=zeros(nx*nx,episode,num_BS);
nablakhistoryBS=zeros(nx*nu,episode,num_BS);
PhistoryBS=zeros(nx*nu,episode,num_BS);
r=0.05;
n=1000;
l=100;
ckBS=zeros(episode,num_BS);
for counter=1:num_BS
    Sigma_w=AmpNoise_forloop_BS(counter)*eye(3);
    kcurrentBS=K_initial;
    Sigmawstar=dlyap((A+B*Kstar),Sigma_w);
    cKstar=trace(Pstar*eye(3));
    for i=1:episode
        disp(['BS ','Num ',num2str(counter), ' Episode ', num2str(i)])
        try
        khistoryBS(:,i,counter)=kcurrentBS(:);
        P=dlyap((A+B*kcurrentBS)',Q+kcurrentBS'*R*kcurrentBS);
        PhistoryBS(:,i,counter)=P(:);
        ckBS(i,counter)=trace(P*eye(3));
        [nablaK,SigmaK,errorgradientBS(i,counter),errorcovarianceBS(i,counter)]=estimationBS(kcurrentBS,l,r,n,A,B,Q,R,Sigma_w,Sigma_0,nx,nu);
        SigmahistoryBS(:,i,counter)=SigmaK(:);
        nablakhistoryBS(:,i,counter)=nablaK(:);
        kcurrentBS=kcurrentBS-varying(counter)*eta_BS(counter)*nablaK;
        catch
            break;
        end
    end
end
% Figure 
BS5=zeros(episode,1);
for i=1:episode
    BS5(i,1)=(ckBS(i,1)-cKstar)/(cKstar);
end
figure;
loglog(BS5,'LineStyle',':','Color','red','LineWidth',1.5)
hold on
ylim([0.01,2])
xlim([0,episode])
grid on
ylabel('$\frac{C{(\hat{K}_i)}-C{(K^*)}}{C{(K^*)}}$','interpreter','latex','FontSize',12)
xlabel('$\mathrm{iteration}~i$','interpreter','latex','FontSize',12)
legend('$\mathrm{GD}~\Sigma_w=10^{-2}I~ \eta=1.08^i$','interpreter','latex','FontSize',9)
save("varying.mat")


%% Required functions
function [nablaK,SigmaK,errorgradient,errorcovariance]=estimationBS(kcurrent,l,r,n,A,B,Q,R,Sigma_w,Sigma_0,nx,nu)
    % true kernal P_K
    P=dlyap((A+B*kcurrent)',Q+kcurrent'*R*kcurrent);
    % true covariance Sigma_K
    Pw=dlyap((A+B*kcurrent),Sigma_w);
    C=2*((R+B'*P*B)*kcurrent+B'*P*A)*Pw;
    ckhistory=zeros(nx,nu);
    xkhistory=zeros(nx,nx);
    for j=1:n
        %disp([' Sample ', num2str(j)])
        xhistory=zeros(nx,l+1);
        xhistory(:,1)=mvnrnd(zeros(1,3),Sigma_0,1);
        X=randommatrix(kcurrent,r);
        K=kcurrent+X;
        ck=0;
        xk=zeros(nx,nx);
        for p=1:l
            xhistory(:,p+1)=realsystem(xhistory(:,p),K*xhistory(:,p),A,B,Sigma_w);
            ck=ck+xhistory(:,p)'*(Q+K'*R*K)*xhistory(:,p);
            xk=xk+xhistory(:,p)*xhistory(:,p)';
        end
        ckhistory=ckhistory+nx*nu/(r^2)*1/l*ck*X;
        xkhistory=xkhistory+1/l*xk;
    end
    nablaK=ckhistory/n;
    errorgradient=norm(nablaK-C);
    SigmaK=xkhistory/n;
    errorcovariance=norm(SigmaK-Pw);
end

% system dynamic
function xnew=realsystem(state,controlinput,A,B,Sigma_w) 
    xnew=A*state+B*controlinput+mvnrnd(zeros(3,1),Sigma_w,1)';
end
% generete the random matrix for exploration
function X=randommatrix(K,r)
    [a,b]=size(K);
    target_norm=r;
    X=randn(a,b);
    current_norm=norm(X,'fro');
    X=X*(target_norm/current_norm);
end
