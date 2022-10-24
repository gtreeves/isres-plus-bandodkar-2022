function J = jac(t,Y,params,~,P,Vn,Vc,An,Am,stage,deval_soln,idx)
% This function calculates the Jacobian of the functions from 
% ftn_dsat with respect to the State Variables


if exist('deval_soln','var')
	Y = deval(deval_soln,t);
	repmat_jac = true;
else
	repmat_jac = false;
end
if exist('idx','var')
	if repmat_jac
		nrep = length(idx);
	end
else
	if repmat_jac
		nrep = length(params);
	end
end


M = length(Y)/6;
I = speye(M);           % Matrices for use later
Z = sparse(M,M);

%
% Unpack state variables
%
un = Y(1:M);
uc = Y(M+1:2*M);
% wn = Y(2*M+1:3*M);
wc = Y(3*M+1:4*M);
vn = Y(4*M+1:5*M);
vc = Y(5*M+1:6*M);


%
% Unpack params
%
params = num2cell(params);
[lambdaU, lambdaW, lambdaV, sigmaU, sigmaW, sigmaV, muU, muW, muV,...
	gamma, psi, alpha, beta, phi, beta0, kappa, delta, K] = params{:};



%
% Make matrix versions of the State Variables
%
% Un  = spdiags(un,0,M,M);
% Uc  = spdiags(uc,0,M,M);
% % Wn = spdiags(wn,0,M,M);
% % Wc = spdiags(wc,0,M,M);
% Vn1 = spdiags(vn,0,M,M);
% Vc1 = spdiags(vc,0,M,M);


spDiag = spdiags(ones(M,1),0,M,M);
Un  = un.*spDiag;
Uc  = uc.*spDiag;
%Wn  = Wn.*spDiag;
%Wc  = Wc.*spDiag;
Vn1 = vn.*spDiag;
Vc1 = vc.*spDiag;



%
% Defining other derivatives
%
% 1. Derivative of the Michaelis-Menten function
% ftn = g*wc./(kappa + wc)
x = linspace(0,1,M)';
g = exp(-x.^2/2/phi^2);
G = g.*kappa./(kappa + wc).^2;
%G = spdiags(G,0,M,M);
G  = G.*spDiag;


% 2. Derivative of the function that produces Cact neg Fbk
% f = ((K*un + 1).^5 - 1)./(K*un + 1).^5;
% f = un.^nH./(K.^nH + un.^nH);
nH = 20;
% F = (5/K*(un/K + 1).^4) ./ (un/K + 1).^10;
F = nH*un.^(nH-1)*K.^nH./(K.^nH + un.^nH).^2;
%F = spdiags(F,0,M,M);
F = F.*spDiag;


% Groupings
if strcmp(stage,'interphase')
	an = An/Vn;
elseif strcmp(stage,'mitosis')
	an = 0;
end




%
% Jacobian elements
% 
% 1st row block, derivatives of f1, the eqn for un.
% f1 = (sigmaU*uc - muU*un) - gamma*un.*vn + beta0*wn;
f1_un = -muU*an*I - gamma*Vn1;
f1_uc = sigmaU*an*I;
f1_wn = beta0*I;
f1_wc = Z;
f1_vn = -gamma*Un;
f1_vc = Z;

J1 = [f1_un f1_uc f1_wn f1_wc f1_vn f1_vc];



% Second row block, derivatives of f2, the eqn for uc.
% f2 = lambdaU*P*uc - (gamma*uc.*vc - beta*g*wc./(kappa + wc) - beta0*wc) ...
%	- (sigmaU*uc - muU*un);
f2_un = muU*An/Vc*I;
f2_uc = lambdaU*Am/Vc*P - gamma*Vc1 - sigmaU*An/Vc*I;
f2_wn = Z;
f2_wc = beta*G + beta0*I;
f2_vn = Z;
f2_vc = -gamma*Uc;

J2 = [f2_un f2_uc f2_wn f2_wc f2_vn f2_vc];




% Third row block, derivatives of f3, the eqn for wn.
% f3 = (sigmaW*wc - muW*wn) + gamma*un.*vn - beta0*wn;
f3_un = gamma*Vn1;
f3_uc = Z;
f3_wn = -muW*an*I - beta0*I;
f3_wc = sigmaW*an*I;
f3_vn = gamma*Un;
f3_vc = Z;

J3 = [f3_un f3_uc f3_wn f3_wc f3_vn f3_vc];




% Fourth row block, derivatives of f4, the eqn for wc.
% f4 = lambdaW*P*wc + (gamma*uc.*vc - beta*g*wc./(kappa + wc) - beta0*wc)...
%	- (sigmaW*wc - muW*wn);
f4_un = Z;
f4_uc = gamma*Vc1;
f4_wn = muW*An/Vc*I;
f4_wc = lambdaW*Am/Vc*P - sigmaW*An/Vc*I - beta*G - beta0*I;
f4_vn = Z;
f4_vc = gamma*Uc;

J4 = [f4_un f4_uc f4_wn f4_wc f4_vn f4_vc];




% Fifth row block, derivatives of f5, the eqn for zn.
% f5 = (sigmaV*vc - muV*vn) - psi*(gamma*un.*zn - beta0*wn);

f5_un = -psi*gamma*Vn1;
f5_uc = Z;
f5_wn = psi*beta0*I;
f5_wc = Z;
f5_vn = -muV*an*I - psi*gamma*Un;
f5_vc = sigmaV*an*I;

J5 = [f5_un f5_uc f5_wn f5_wc f5_vn f5_vc];




% Sixth row block, derivatives of f6, the eqn for zc.
% f6 = lambdaV*P*vc - psi*(gamma*uc.*vc - beta*g*wc./(kappa + wc)) + ...
%	1 + delta*f - alpha*vc - (sigmaV*vc - muV*vn);
f6_un = delta*F;
f6_uc = -psi*gamma*Vc1;
f6_wn = Z;
f6_wc = psi*(beta*G + beta0*I);
f6_vn = muV*An/Vc*I;
f6_vc = lambdaV*Am/Vc*P - psi*gamma*Uc - alpha*I - sigmaV*An/Vc*I;

J6 = [f6_un f6_uc f6_wn f6_wc f6_vn f6_vc];

J = [J1; J2; J3; J4; J5; J6];


%
% Do we repeat the jacobian?
%
if repmat_jac
	C = repmat({J},1,nrep);
	J = blkdiag(C{:});	
end


