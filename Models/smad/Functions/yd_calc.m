function ycalc_d = yd_calc(k,y0)
load('tgf.mat');
TGF1 = transpose(tgf); 
tspan1 = TGF1(1,1:15);

ycalc_d = RtoODE1(k,tspan1,y0);

end
