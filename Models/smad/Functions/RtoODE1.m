function solpts = RtoODE1(k,tspan,y0)
solpts = RtoODE(k,tspan,y0);
solpts = sum([solpts(17,:);solpts(19,:);solpts(22,:);solpts(24,:);2*solpts(25,:)]); % EGFP Smad2 in the nucleus

end