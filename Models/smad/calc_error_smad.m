function [F,Phi] = error_calc_smad(k,model)

k     = 10.^(k); %not sure why
nsets = size(k,1);
F     = zeros(nsets,1);
Phi   = zeros(nsets,1); %no penalty
% initial conditions
y0 = [1 0.066 0 0 10000 60.6 60.6 0 0 50.8 0 0 0 0 0 28.5 28.5 0 0 50.8 0 0 0 0 0];

load('error_bar_1.mat')
load('error_bar_2.mat')
errorbar = transpose(errorbar1.error_upper - errorbar2.error_lower);
errorbar_left = errorbar(1,1:15);
errorbar_right = errorbar(1,16:35);
% figure 2D (before SB addition)
load('smad.mat');
SMAD2 = transpose(smad);
yvalsd = SMAD2(1,1:15);
load('tgf.mat');
TGF1 = transpose(tgf); 
tspan1 = TGF1(1,1:15);
% figure 2D (after SB addition)
yvalsd_SB = SMAD2(1,16:35);



parfor i=1:nsets
    ktgf = k(i,1);
    kphos = k(i,2);
    konSB = k(i,3);
    kon = k(i,4);
    kdephos = k(i,5);
    CIF = k(i,6);
    kin = k(i,7);
    kex = k(i,8);
    koff = k(i,9);
    koffSB = k(i,10);
    
    ka =[ktgf kphos konSB kon kdephos CIF kin kex koff koffSB];
    

    % calculated value of figure 2D (before SB addition)
    ycalc_d = yd_calc(ka,y0);
    % calculated value of figure 2D (before SB addition)
    inc = RtoODE(ka,tspan1,y0);
    ySB = inc(:,15);
    ycalc_d_SB = yd_calcSB(ka,ySB);


    F(i) = sum(((ycalc_d - yvalsd)./errorbar_left).^2,'all')+sum(((ycalc_d_SB - yvalsd_SB)./errorbar_right).^2,'all'); %+ sum((ycalc_e(:,:,i) - yvalse).^2) + sum((ycalc_f(:,:,i) - yvalsf).^2);
    
end

end 
