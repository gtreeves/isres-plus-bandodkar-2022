function dVdt = ftn_manueqns(t,V,params)

%unpackstruct(params);

% Unpack state variables
M   = length(V)/4;
Hb  = V(1:M);
Kr  = V(M+1:2*M);
Gt  = V(2*M+1:3*M);
Kni = V(3*M+1:4*M);


% Make x
x  = linspace(params.xL,params.xU,M)';
dx = x(2) - x(1);


% Make Bcd, Cad, Tll
Bcd  = params.A*exp(-params.phi*x);
Cad1 = interp1(params.tc,params.Cad',t)';
Cad  = interp1(params.xmat,Cad1,x);
Tll1 = interp1(params.tc,params.Tll',t)';
Tll  = interp1(params.xmat,Tll1,x);
Hkb1 = interp1(params.tc,params.Hkb',t)';
Hkb  = interp1(params.xmat,Hkb1,x);



%
% make P-matrix
%
P       = ones(M,1)*[1 -2 1]; 
P       = spdiags(P,[-1 0 1],M,M);
P(1,2)  = 2*P(1,2);     % reflective BC at the RHS    
% absorbing BC at i = N means no modification there.
if params.xU == 1
	P(end,end-1) = 2*P(end,end-1);
end

if t > params.t13m && t < params.t14i
	params.R = 0*params.R;
end
if M == params.M13
	x14 = linspace(params.xL,params.xU,params.M14)';
	dx14 = x(2) - x(1);
	params.D = params.D*dx14^2/dx^2;
end


%
% Regulatory variables
%
V2  = [Bcd Cad Tll Hb Kr Gt Kni Hkb];
uHb  = V2*params.T(1,:)' + params.ha;
uKr  = V2*params.T(2,:)' + params.ha;
uGt  = V2*params.T(3,:)' + params.ha;
uKni = V2*params.T(4,:)' + params.ha;


%
% diffyQ's
%
dHbdt  = params.R(1)*gftn(uHb)  + params.D(1)*P*Hb  - params.lambda(1)*Hb;
dKrdt  = params.R(2)*gftn(uKr)  + params.D(2)*P*Kr  - params.lambda(2)*Kr;
dGtdt  = params.R(3)*gftn(uGt)  + params.D(3)*P*Gt  - params.lambda(3)*Gt;
dKnidt = params.R(4)*gftn(uKni) + params.D(4)*P*Kni - params.lambda(4)*Kni;


%
% Package dVdt
%
dVdt = [dHbdt; dKrdt; dGtdt; dKnidt];

end


function g = gftn(u)
    g = 0.5*(u./sqrt(u.^2+1) + 1);
end


