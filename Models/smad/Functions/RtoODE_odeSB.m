function solpts = RtoODE_odeSB(k,tspan,y0)
sol = ode15s(@(t,y)diffunSB(t,y,k),tspan,y0); %45
solpts = deval(sol,tspan);
end