% ACFD reading and plotting data from c++ code
Y(1)=0;Y(102)=1;X(1)=0;X(102)=1;
Y(2:101)=0.005:0.01:0.995;
X(2:101)=0.005:0.01:0.995;
XX=X(2:101);YY=Y(2:101);
filename = 'u_v_data.txt';
delimiterIn = ',';
headerlinesIn = 0;
uvdata_quick = importdata(filename,delimiterIn,headerlinesIn);
filename = 'u_data.txt';
delimiterIn = ',';
headerlinesIn = 0;
udata_quick = importdata(filename,delimiterIn,headerlinesIn);
filename = 'v_data.txt';
delimiterIn = ',';
headerlinesIn = 0;
vdata_quick = importdata(filename,delimiterIn,headerlinesIn);
filename = 'temperature.txt';
delimiterIn = ',';
headerlinesIn = 0;
temperature_quick = importdata(filename,delimiterIn,headerlinesIn);

figure(1)
plot(uvdata_upwind(1,:),Y,'k','Linewidth',1.5);
hold on
plot(uvdata_quick(1,:),Y,'r','Linewidth',1.5);
% plot(URE100CHOI(:,1),URE100CHOI(:,2),'o','Markeredgecolor','blue');
% plot(UGhiaetalRE100(:,1),UGhiaetalRE100(:,2),'s','Markeredgecolor','g','Markersize',8);
ylabel('Y/L','Fontweight' ,'bold');
xlabel('U/U*','Fontweight' ,'bold');
title({'U velocity variation along center vertical line'});
legend('upwind','quick','Location','northwest');
grid on
% ylim([0 1.01]);
hold off

figure(2)
plot(X,uvdata_upwind(2,:),'k','Linewidth',1.5);
hold on
plot(X,uvdata_quick(2,:),'r','Linewidth',1.5);
% plot(VchoietalRE100(:,1),VchoietalRE100(:,2),'o','Markeredgecolor','blue');
% plot(VRE100ghia(:,1),VRE100ghia(:,2),'s','Markeredgecolor','g','Markersize',8);
xlabel('X/L','Fontweight' ,'bold');
ylabel('V/U*','Fontweight' ,'bold');
title({'V velocity variation along center horizontal line'});
legend('upwind','quick');
% xlim([-0.01 1.01]);
grid on
hold off


figure(3)
contourf(X,Y,udata_upwind,20);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12} U velocity -- UPWIND'});
colorbar;

figure(4)
contourf(X,Y,vdata_upwind,20);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12} V velocity -- UPWIND'});
colorbar;

figure(5)
contourf(X,Y,temperature_upwind,20);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12} Temperature - Upwind'});
colorbar;


figure(6)
contourf(X,Y,udata_quick,20);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12} U velocity -- QUICK'});
colorbar;

figure(7)
contourf(X,Y,vdata_quick,20);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12} V velocity -- QUICK'});
colorbar;


figure(8)
contourf(X,Y,temperature_quick,20);
xlabel('X/L','Fontweight' ,'bold','Fontsize',12);
ylabel('Y/L','Fontweight' ,'bold','Fontsize',12);
title({'\fontsize{12} Temperature - Quick'});
colorbar;

Nul_upwind = -(temperature_upwind(2:101,2)-temperature_upwind(2:101,1))/0.005;
Nur_upwind = -(temperature_upwind(2:101,102)-temperature_upwind(2:101,101))/0.005;
Nul_quick = -(temperature_quick(2:101,2)-temperature_quick(2:101,1))/0.005;
Nur_quick = -(temperature_quick(2:101,102)-temperature_quick(2:101,101))/0.005;
figure(9)
plot(Nul_upwind,YY,'k','Linewidth',1.5);
hold on
plot(Nul_quick,YY,'g','Linewidth',1.5);
plot(Nu_hot(:,1),Nu_hot(:,2),'o','Markeredgecolor','blue');
ylabel('Y/L','Fontweight' ,'bold');
xlabel('Nu','Fontweight' ,'bold');
title({'Nusselt number variation along X=0 line (Hot wall)'});
legend('Upwind','Quick','Benchmark');
grid on
hold off


figure(10)
plot(Nur_upwind,YY,'r','Linewidth',1.5);
hold on
plot(Nur_quick,YY,'b','Linewidth',1.5);
plot(Nu_cold(:,1),Nu_cold(:,2),'s','Markeredgecolor','g','Markersize',8);
ylabel('Y/L','Fontweight' ,'bold');
xlabel('Nu','Fontweight' ,'bold');
title({'Nusselt number variation along X=1 line (Cold wall)'});
legend('Upwind','quick','Benchmark');
grid on
hold off