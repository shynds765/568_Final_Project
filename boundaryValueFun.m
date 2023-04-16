function Psi = boundaryValueFun(lamda0,Tf,X0,mu,debris_sol,ref_sol)
    options_ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);
    Xug0 = [X0;lamda0];

    Q = diag(eye(6));

    R = diag(eye(3));

    d0 = 10; % km

    % System Parameters 
    gamma = 0;
    rho = 1;
    [~,Xaug] = ode45(@(t,Xaug) augEOM(t,Xaug,Q,R,debris_sol,ref_sol,mu,d0,gamma,rho),[0,Tf], Xug0, options_ode);
    lambda_f = Xaug(end,7:12);

    % Calc Boundary Condition
    % Costates have to be 0

    % Psi = 0
    Psi = lambda_f;
end