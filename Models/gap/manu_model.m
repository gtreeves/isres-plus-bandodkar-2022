function [expNC13, expNC14]=manu_model(p)

load('Models/gap/Mats/gapgrn.mat');   % loads a variable called params

params.R        = p(1:4);
params.D        = p(5:8);
params.lambda   = p(9:12);
params.ha       = -2.5;
params.T(1,:)   = p(13:20);
params.T(2,:)   = p(21:28);
params.T(3,:)   = p(29:36);
params.T(4,:)   = p(37:44);

M13 = params.M13; %raz
M14 = params.M14; %raz
xL  = params.xL;  %raz
xU  = params.xU;  %raz


%
% Setup initial conditions for each gene for nc13
%
x13  = linspace(xL,xU,M13)';   %linspace(params.xL,params.xU,M13)'; %raz
x1   = params.Hb_data.x;
y1   = params.Hb_data.y;
Hb0  = interp1(x1,y1,x13);      Hb0(isnan(Hb0)) = 0; % setting nans to zero.
Kr0  = zeros(M13,1); 
Gt0  = Kr0; 
Kni0 = Kr0;
V130 = [Hb0; Kr0; Gt0; Kni0];


%
% Run nc 13
%
ftnhandle = @ftn_manueqns;
tspan13   = [0, params.t14i];
[T13,V13] = ode15s(ftnhandle,tspan13,V130,[],params);


%
% Extract gene expression information for nc13
%
V13     = 0.8*V13;
Hb13    = V13(:,1:M13);
Kr13    = V13(:,M13+1:2*M13);
Gt13    = V13(:,2*M13+1:3*M13);
Kni13   = V13(:,3*M13+1:4*M13);
expNC13 = [Hb13(end,:)', Kr13(end,:)', Gt13(end,:)', Kni13(end,:)'];


%
% Setup initial conditions for each gene for nc14
%
x14     = linspace(xL,xU,M14)';
Hb140   = interp1(x13,Hb13(end,:)',x14);
Kr140   = interp1(x13,Kr13(end,:)',x14);
Gt140   = interp1(x13,Gt13(end,:)',x14);
Kni140  = interp1(x13,Kni13(end,:)',x14);
V140    = [Hb140; Kr140; Gt140; Kni140];


%
% Run nc 14
%
tspan14   = [params.t14i, params.tc(3:end)];
[T14,V14] = ode15s(ftnhandle,tspan14,V140,[],params);


%
% Extract gene expression information for nc14
%
Hb14    = V14(:,1:M14);
Kr14    = V14(:,M14+1:2*M14);
Gt14    = V14(:,2*M14+1:3*M14);
Kni14   = V14(:,3*M14+1:4*M14);

% We do not have all Time Class data for Kr and Kni;
expNC14 = [Hb14(2:end,:)', Kr14([2 3 5 7 8 9],:)', Gt14(2:end,:)', Kni14([2 3 4 6 8 9],:)'];

end