function fout = colebrook_solve(x)
    epsilon = x(1);
    Dh = x(2);
    rho = x(3);
    velocity = x(4);
    mu = x(5);
    Re = Dh * rho * velocity/mu;
    %Haaland equation as initial guess.
    haaland_LHS = -1.8*log10((epsilon/Dh/3.7)^1.11 + 6.9/Re);
    x0 = (1/haaland_LHS)^2;
    fout = fzero(@colebrook_zero,x0);
    
    function deltaout = colebrook_zero(f)
        deltaout = 1/sqrt(f) + 2*log10(epsilon/(3.7*Dh) + 2.51/(Re*sqrt(f)));
    end
end
