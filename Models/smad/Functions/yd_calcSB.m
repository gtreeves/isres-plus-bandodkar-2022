function ycalc_d_SB = yd_calcSB(k,y0)
load('tgf.mat');
TGF1 = transpose(tgf); 
tspanSB = TGF1(1,16:35);

ycalc_d_SB = RtoODESB(k,tspanSB,y0);

end
