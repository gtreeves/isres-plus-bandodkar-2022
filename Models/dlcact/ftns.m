 function dYdt = ftns(~,Y,params,x,P,Vn,Vc,An,Am,stage)
% This is the function-handle where all the differential equations reside

M = length(Y)/6;

%
% Unpack state variables
%
un = Y(1:M);
uc = Y(1*M+1 :2*M);
wn = Y(2*M+1 :3*M);
wc = Y(3*M+1 :4*M);
vn = Y(4*M+1 :5*M);
vc = Y(5*M+1 :6*M);

%
% Unpack params
%
params = num2cell(params);
[lambdaU, lambdaW, lambdaV, sigmaU, sigmaW, sigmaV, muU, muW, muV,...
	gamma, psi, alpha, beta, phi, beta0, kappa, delta, K] = params{:};




%
% Defining Toll and exchange-related terms
%
g   = exp(-0.5*(x/phi).^2);
BEn = gamma*un.*vn - beta0*wn;                                % dl/Cact binding eqm in nucleus
BEc = gamma*uc.*vc - beta*g.*wc./(kappa + wc) - beta0*wc;     % dl/Cact binding eqm in cytoplasm
%BEc = gamma*uc.*vc - beta*g.*wc - beta0*wc;       
if strcmp(stage,'interphase')
% 	f = ((un/K + 1).^5 - 1)./(un/K + 1).^5;
	nH = 20;
	f = un.^nH./(K.^nH + un.^nH);
	nucexch_u = An*(sigmaU *uc- muU *un);  % dl nuc exchange
	nucexch_w = An*(sigmaW *wc- muW *wn);  % dl/Cact nuc exchange
	nucexch_v = An*(sigmaV *vc- muV *vn);  % Cact nuc exchange
    
elseif strcmp(stage,'mitosis')
	f = 0;
	nucexch_u = 0;
	nucexch_w = 0;
	nucexch_v = 0;
end



%
% Define differential equations
%
f1 = nucexch_u/Vn - BEn;
f2 = lambdaU*Am/Vc*P*uc - nucexch_u/Vc - BEc;
f3 = nucexch_w/Vn + BEn;
f4 = lambdaW*Am/Vc*P*wc - nucexch_w/Vc + BEc;
f5 = nucexch_v/Vn - psi*BEn;
f6 = lambdaV*Am/Vc*P*vc - nucexch_v/Vc - psi*BEc + 1/Vc + delta*f - alpha*vc;

if strcmp(stage,'mitosis')
	f1 = zeros(M,1);
	f3 = zeros(M,1);
	f5 = zeros(M,1);
end

dYdt = [f1; f2; f3; f4; f5; f6];


if any(isnan(dYdt(:))) || any(isinf(dYdt(:))) || any(isnan(Y(:)))
	1;
end



