function solpts = RtoODE(k,tspan,y0)
sol = ode15s(@(t,y)diffun(t,y,k),tspan,y0); %45
solpts = deval(sol,tspan);
end