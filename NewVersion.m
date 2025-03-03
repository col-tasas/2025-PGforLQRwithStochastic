%% Parameters 
% define the true system and the weighting matrix 
A=[1.01,0.01,0;0.01,1.01,0.01;0,0.01,1.01];
B=1*eye(3); % input matrix B
nx=size(A,2); % number of states
nu=size(B,2); % number of control inputs
Q=0.001*eye(3); % weighting matrix Q
R=1*eye(3); % weighting matrix R
Sigma_w=1*eye(3);
Sigma_0=0.001*eye(3);
% K_initial=diag([-1.5,-1,-0.5]);
[Pstar1,K_initial]=idare(A,B,50*Q,R);
K_initial=-K_initial;
[Pstar,Kstar]=idare(A,B,Q,R);
Kstar=-Kstar;
CKstar=trace(Pstar*Sigma_w);
iteration=100;
number=10000;
CKPG=zeros(number,iteration);
CKPG2=zeros(number,iteration);
CKPG3=zeros(number,iteration);
CKPGV=zeros(1,iteration);
KPG=zeros(size(K_initial,1),size(K_initial,2),number);
KPG2=zeros(size(K_initial,1),size(K_initial,2),number);
KPG3=zeros(size(K_initial,1),size(K_initial,2),number);
KPGnoisefree=K_initial;
for i=1:number
    KPG(:,:,i)=K_initial;
    KPG2(:,:,i)=K_initial;
    KPG3(:,:,i)=K_initial;
end
noisevariancce=0.6;
for i=1:iteration
    i
    stepsie=0.12;
     for j=1:number
    P=dlyap((A+B*KPG(:,:,j))',Q+KPG(:,:,j)'*R*KPG(:,:,j));
    Pw=dlyap((A+B*KPG(:,:,j)),Sigma_w);
    CKPG(j,i)=trace(P*Sigma_w);
    %PGstep=hck(CKPG(1,1),B,Q,R,Sigma_w,CKstar);
    
    KPG(:,:,j)=KPG(:,:,j)-stepsie*(((R+B'*P*B)*KPG(:,:,j)+B'*P*A)*Pw+normrnd(0,noisevariancce,[3,3]));
     end
    % PGV
    P=dlyap((A+B*KPGnoisefree)',Q+KPGnoisefree'*R*KPGnoisefree);
    Pw=dlyap((A+B*KPGnoisefree),Sigma_w);
    CKPGV(1,i)=trace(P*Sigma_w);
    %PGVstep=hck(CKPG(1,i),B,Q,R,Sigma_w,CKstar);
    KPGnoisefree=KPGnoisefree-stepsie*((R+B'*P*B)*KPGnoisefree+B'*P*A)*Pw;
     for j=1:number
    P=dlyap((A+B*KPG2(:,:,j))',Q+KPG2(:,:,j)'*R*KPG2(:,:,j));
    Pw=dlyap((A+B*KPG2(:,:,j)),Sigma_w);
    CKPG2(j,i)=trace(P*Sigma_w);
    %PGstep=hck(CKPG(1,1),B,Q,R,Sigma_w,CKstar);
    stepsie=0.01;
    KPG2(:,:,j)=KPG2(:,:,j)-stepsie*(((R+B'*P*B)*KPG2(:,:,j)+B'*P*A)*Pw+normrnd(0,noisevariancce,[3,3]));
     end
     stepsie=0.12;
     for j=1:number
    P=dlyap((A+B*KPG3(:,:,j))',Q+KPG3(:,:,j)'*R*KPG3(:,:,j));
    Pw=dlyap((A+B*KPG3(:,:,j)),Sigma_w);
    CKPG3(j,i)=trace(P*Sigma_w);
    %PGstep=hck(CKPG(1,1),B,Q,R,Sigma_w,CKstar);
    KPG3(:,:,j)=KPG3(:,:,j)-stepsie*(((R+B'*P*B)*KPG3(:,:,j)+B'*P*A)*Pw+normrnd(0,0.03,[3,3]));
     end

end
figure;
CKstarfigure=ones(1,iteration)*CKstar;
%loglog(1:1:iteration,(CKPI-CKstarfigure)/CKstar,'LineStyle','-','Color','red','LineWidth',1.5);
%hold on
%loglog(1:1:iteration,(CKNPGV-CKstarfigure)/CKstar,'LineStyle','-','Color','blue','LineWidth',1.5)
%hold on
%loglog(1:1:iteration,(CKNPG-CKstarfigure)/CKstar,'LineStyle','-','Color','m','LineWidth',1.5)
%hold on
loglog(1:1:iteration,(CKPGV-CKstarfigure)/CKstar,'LineStyle',':','Color','black','LineWidth',1.5)
hold on
loglog(1:1:iteration,(abs(mean(CKPG3,1)-CKstarfigure))/CKstar,'LineStyle','--','Color','blue','LineWidth',1.5)
hold on
loglog(1:1:iteration,(abs(mean(CKPG,1)-CKstarfigure))/CKstar,'LineStyle','-.','Color','green','LineWidth',1.5)
hold on
loglog(1:1:iteration,(abs(mean(CKPG2,1)-CKstarfigure))/CKstar,'LineStyle','-','Color','red','LineWidth',1.5)
hold on

ylim([0.01,50])
xlim([0,100])
grid on
ylabel('$\frac{C(\hat{K}_i)-C(K^*)}{C(K^*)}$','interpreter','latex','FontSize',12)
xlabel('$\mathrm{iteration}~i$','interpreter','latex','FontSize',12)
legend('$\sigma=0~\eta=0.12$','$\sigma=0.03~\eta=0.12$','$\sigma=0.6~\eta=0.12$','$\sigma=0.6~\eta=0.01$','interpreter','latex','FontSize',9)
