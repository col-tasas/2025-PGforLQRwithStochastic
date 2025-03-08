load("NPGVarying01261test2.mat")
BS1=zeros(100,1);
BS2=zeros(100,1);
BS3=zeros(100,1);
BS4=zeros(100,1);
for i=1:100
    % BS1(i,1)=(ckBS(i,5)-cKstar)/(cKstar);
    % BS2(i,1)=(ckBS(i,8)-cKstar)/(cKstar);
    % BS3(i,1)=(ckBS(i,12)-cKstar)/(cKstar);
    % BS4(i,1)=(ckFA(i,3)-cKstar)/cKstar;
    BS1(i,1)=(ckBS(i,9)-cKstar)/(cKstar);
    BS2(i,1)=(ckBS(i,12)-cKstar)/(cKstar);
    BS3(i,1)=(ckBS(i,28)-cKstar)/(cKstar);
    %BS4(i,1)=(ckFA(i,3)-cKstar)/cKstar;
end
load("NPGVaryingtest1920.mat")
for i=1:100
    BS4(i,1)=(ckFA(i,1)-cKstar)/cKstar;
end

figure;
%plot(FA1,'LineStyle','-','Color','blue','LineWidth',1.5)
%hold on
%plot(FA2,'LineStyle','--','Color','red','LineWidth',1.5)
%hold on
%plot(FA3,'LineStyle',':','Color','black','LineWidth',1.5)
%hold on

loglog(BS1,'LineStyle',':','Color','red','LineWidth',1.5)
hold on
loglog(BS2,'LineStyle','-','Color','green','LineWidth',1.5)
hold on
loglog(BS3,'LineStyle','-.','Color','blue','LineWidth',1.5)
hold on
loglog(BS4,'LineStyle','--','Color','black','LineWidth',1.5)
hold on

ylim([0.01,2])
xlim([0,100])
grid on
ylabel('$\frac{C{(\hat{K}_i)}-C{(K^*)}}{C{(K^*)}}$','interpreter','latex','FontSize',12)
xlabel('$\mathrm{iteration}~i$','interpreter','latex','FontSize',12)
legend('$\mathrm{NPG}~\Sigma_w=10^{-4}I~ \eta=0.1$','$\mathrm{NPG}~\Sigma_w=10^{-2}I~ \eta=0.1$','$\mathrm{NPG}~\Sigma_w=10^{0}I~ \eta=0.1$','$\mathrm{NPG}~\Sigma_w=10^{-2}I~ \eta_i=1.01^i\eta$','interpreter','latex','FontSize',9)