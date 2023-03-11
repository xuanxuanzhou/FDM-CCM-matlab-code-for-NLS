%stab  chebyshevstab
%专门画p<=5图(plot the figures for p<=5)
tau=0.001;
Tt=10;
Nt=1/tau*Tt;
tp=linspace(0,Tt,Nt+1);
%第一组图 
figure;
% axes('Position',[0.13,0.13,0.8,0.8]);

plot(tp,log10(real(unum(1,:))),'r-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp,log10(real(unum(2,:)) ),'k-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp,log10(real(unum(3,:)) ),'b-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp,log10( real(unum(4,:)) ),'g-.','MarkerSize',12,'LineWidth',1.5)
% hold on
% plot(tp,real(unum(5,:)),'k.','MarkerSize',12,'LineWidth',1.5)
% hold on
% plot(tp,real(unum(6,:)),'r.','MarkerSize',12,'LineWidth',1)
% hold on
% plot(tp,real(unum(7,:)),'b.','MarkerSize',12,'LineWidth',1)
xlabel('t','fontsize',16);ylabel('log10(e(\omega,\eta))','fontsize',16);
title('p=5,\omega=\lambda^2/4+10','fontsize',16);
legend('\eta=4*10^{-1}','\eta=4*10^{-2}','\eta=4*10^{-3}','\eta=4*10^{-4}');


for k=1:100:Nt+1
    tp1((k-1)/100+1,1)=tp(k);
    ustab21(1,(k-1)/100+1,1)=ustab2(1,k);
    ustab21(2,(k-1)/100+1,1)=ustab2(2,k);
    ustab21(3,(k-1)/100+1,1)=ustab2(3,k);
    ustab21(4,(k-1)/100+1,1)=ustab2(4,k);
end
%% 以下两行嵌入局部放大图用
%hold on
%axes('Position',[0.3,0.25,0.44,0.45]);
%% 以下图太密看不清
figure;
plot(tp,log10( real(ustab2(1,:)) ),'r-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp,log10( real(ustab2(2,:)) ),'k-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp,log10( real(ustab2(3,:)) ),'b-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp,log10( real(ustab2(4,:)) ),'g-.','MarkerSize',12,'LineWidth',1.5)
figure;
plot(tp1,log10( real(ustab21(1,:)) ),'r-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp1,log10( real(ustab21(2,:)) ),'k-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp1,log10( real(ustab21(3,:)) ),'b-.','MarkerSize',12,'LineWidth',1.5)
hold on
plot(tp1,log10( real(ustab21(4,:)) ),'g-.','MarkerSize',12,'LineWidth',1.5)
xlabel('t','fontsize',16);ylabel('log10(e_{\theta}(\omega,\eta))','fontsize',16);
%axis([0,Tt,0,20])
title('p=5,\omega=\lambda^2/4+10','fontsize',16);
legend('\eta=4*10^{-1}','\eta=4*10^{-2}','\eta=4*10^{-3}','\eta=4*10^{-4}');




