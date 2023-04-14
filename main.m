clear

options_fsolve = optimset('TolFun',1e-8,'TolX',1e-8,'Display','off');
options_ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);

lambda0_guess = ;
Tf = ;
lambda0 = fsolve(@(lamda0) boundaryValueFun(lamda0,Tf,X0),lambda0_guess,options_fsolve);