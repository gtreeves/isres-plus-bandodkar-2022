function [Dl,sc] = dosage_sensitivity(fcn,xb)


ntotal      = 10;
dl0         = 1;
cact0       = 1;
reltol      = 1e-7;
abstol      = 1e-7;
[soln,K]    = fcn(xb,dl0,cact0,reltol,abstol);
xb          = [xb,K];
dlnucwt     = soln.dlNuc.NC14(:,end);
% plot(dlnucwt)
% hold on;

Dl = zeros(51,ntotal+1);
Dl(:,1) = dlnucwt;

sc = zeros(51,ntotal);
for i=1:ntotal
    delx = randrange(1e-4,1e-5,1)*dl0;
    soln = run_model(xb,dl0 + delx,cact0,reltol,abstol);
    dl   = soln.dlNuc.NC14(:,end);
    
    dely = dl - dlnucwt;
    delx = delx - 0;
    
    sc(:,i) = dely./dlnucwt/delx/dl0;
%     plot(dl)
    Dl(:,i+1) = dl;
end



end

