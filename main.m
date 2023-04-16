clear

options_fsolve = optimset('TolFun',1e-8,'TolX',1e-8,'Display','off');
options_ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);



lambda0_guess = zeros(6,1);
Tf = ;

a0 = 6.7862E3; % km
e0 = 0;
i0 = 0;
RAAN0 = 0;
w0 = 0;
M0 = 0;

X0 = [a0;e0;i0;RAAN0;w0;M0];
lambda0 = fsolve(@(lamda0) boundaryValueFun(lamda0,Tf,X0),lambda0_guess,options_fsolve);