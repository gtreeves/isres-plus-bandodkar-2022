clear
clc
close all

matfolder = './';
files = extractFileLocations(matfolder,'mat');

% use database functions to get files


load(file,'xb')

%
% run manu model
%
[NC13, NC14] = manu_model(xb);

% load exp data
exp = load('Manu/Mats/ExpDataManu.mat');


% load other necessary matfiles 
load('Manu/Mats/gapgrn.mat');   % loads a variable called params
struct2vars(params)

x13 = linspace(xL,xU,M13);
x14 = linspace(xL,xU,M14);



%
% NC13
% Make plots for nc13. We only have data from the final timepoint
%
figure
t = tiledlayout(2,2);
plot(x13,exp.NC13(:,1)); hold on
plot(x13,NC13(:,1))
legend('Exp','Model')
title('Hb')

plot(x13,exp.NC13(:,2)); hold on
plot(x13,NC13(:,2))
legend('Exp','Model')
title('Kr')

plot(x13,exp.NC13(:,3)); hold on
plot(x13,NC13(:,3))
legend('Exp','Model')
title('Gt')

plot(x13,exp.NC13(:,4)); hold on
plot(x13,NC13(:,4))
legend('Exp','Model')
title('Kni')

title(t,'Model vs Exp in NC13')


%
% NC14
% Make plots for nc13. We only have data from the final timepoint
%
thb = 8; tkr = 14; tgt = 22; tkni = 28; 

figure
t = tiledlayout(2,2);
plot(x14,exp.NC14(:,thb)); hold on
plot(x13,NC13(:,thb))
legend('Exp','Model')
title('Hb')

plot(x13,exp.NC13(:,tkr)); hold on
plot(x13,NC13(:,tkr))
legend('Exp','Model')
title('Kr')

plot(x13,exp.NC13(:,tgt)); hold on
plot(x13,NC13(:,tgt))
legend('Exp','Model')
title('Gt')

plot(x13,exp.NC13(:,tkni)); hold on
plot(x13,NC13(:,tkni))
legend('Exp','Model')
title('Kni')

title(t,'Model vs Exp in NC13')



































