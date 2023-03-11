%stab  chebyshevstab
%专门画p>5图,每个p 单独画（plot the figures for p>5）

%% T=100
tau=0.001;
Tt=100;
Nt=1/tau*Tt;
tp=linspace(0,Tt,Nt+1);
%第一组图 
figure;
% axes('Position',[0.13,0.13,0.8,0.8]);

plot(tp, real(unum(1,:)) ,'r-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp, real(unum(2,:)) ,'k-.','MarkerSize',12,'LineWidth',1.5)
xlabel('t','fontsize',16);ylabel('e(\omega,\eta)','fontsize',16);
title('p=6,\eta=10^{-7},\xi_c(6)=0.2923','fontsize',16);
legend('\xi=0.2922','\xi=0.2924');
% title('p=7,\eta=0,\xi_c(7)=0.4409','fontsize',16);
% legend('\xi=0.4407','\xi=0.4411');
% title('p=8,\eta=0,\xi_c(8)=0.5668','fontsize',16);
% legend('\xi=0.5666','\xi=0.5670');

%T=10,100, T=100,500
for k=1:500:Nt+1
    tp1((k-1)/500+1,1)=tp(k);
    ustab21(1,(k-1)/500+1,1)=ustab2(1,k);
    ustab21(2,(k-1)/500+1,1)=ustab2(2,k);
end


figure;
plot(tp1, real(ustab21(1,:)) ,'r-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp1, real(ustab21(2,:)) ,'k-.','MarkerSize',12,'LineWidth',1.5)

xlabel('t','fontsize',16);ylabel('e_{\theta}(\omega,\eta)','fontsize',16);
%axis([0,Tt,0,20])10^{-4}
title('p=6,\eta=10^{-7},\xi_c(6)=0.2923','fontsize',16);
legend('\xi=0.2922','\xi=0.2924');
% title('p=7,\eta=0,\xi_c(7)=0.4409','fontsize',16);
% legend('\xi=0.4407','\xi=0.4411');
% title('p=8,\eta=0,\xi_c(8)=0.5668','fontsize',16);
% legend('\xi=0.5666','\xi=0.5670');


