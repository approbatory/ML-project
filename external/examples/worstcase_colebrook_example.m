function [x_LB,x_UB,f_LB,f_MID,f_UB,minus_percent,plus_percent] = worstcase_colebrook_example
    epsilon = 1e-3; %Relative roughness; should come from Moody chart but meah.
    Dh = 0.05; %Pipe diameter; we'll use 5 cm.
    rho = 1000; %Density of water; kg/m^3
    velocity = 30; %We'll use 30 m/s.
    mu = 0.001; %Water, so we'll use 1 Poise.

    x = [epsilon Dh rho velocity mu]';
    
    %uncertainties = [0.05*1e-3 0.001 0.01*1000 5 1e-4]';
    uncertainties = 0.25*x; %25% error in all values.
    [x_LB,x_UB,f_LB,f_MID,f_UB,minus_percent,plus_percent] = worstcase(@colebrook_solve,x,uncertainties);
end

