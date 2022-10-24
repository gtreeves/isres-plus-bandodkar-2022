%load errorbars
load('error_bar_1.mat')
load('error_bar_2.mat')
errorbar1 = transpose(errorbar1.error_upper);
errorbar2 = transpose(errorbar2.error_lower);
%load experimental data
load('smad.mat');
SMAD2 = transpose(smad);
load('tgf.mat');
TGF1 = transpose(tgf); 
tspan1 = TGF1(1,1:15); tspan2 = TGF1(1,16:35);
time1 = 0:46; time2 = 46:150;
%error bars
pos = errorbar1 - SMAD2;
neg = SMAD2 - errorbar2;
%initial conditions
y0 = [1 0.066 0 0 10000 60.6 60.6 0 0 50.8 0 0 0 0 0 28.5 28.5 0 0 50.8 0 0 0 0 0];

figure
%load a parameter set
k = xb;
Smad2n = RtoODE1(k(1:10),time1,y0);
SBinc = RtoODE(k(1:10),time1,y0);
ySB = SBinc(:,end);
Smad2nSB = RtoODESB(k(1:10),time2,ySB);
y_afterSB = RtoODE_odeSB(k(1:10),time2,ySB);

%plot
errorbar(tspan1,SMAD2(1,1:15),pos(1,1:15),neg(1,1:15),'ro')
hold on
plot(time1,Smad2n,'b-','Linewidth',2)
plot(time2,Smad2nSB,'b-','Linewidth',2)
errorbar(tspan2,SMAD2(1,16:35),pos(1,16:35),neg(1,16:35),'ro')
hold off

legend('Data','ISRES+','','')
ylabel('Nuclear EGFP-Smad2 [nM]')
xlabel('Time [min]')
title('Smad signaling model')
set(gca,'FontSize',14)