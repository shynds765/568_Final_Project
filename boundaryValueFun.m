function Psi = boundaryValueFun(lamda0,Tf,X0,mu,debris_sol)
    options_ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);
    Xug0 = [X0;lamda0];

    Q = eye(6);
    Q(6,6) = 0; % Cost of M is 0 (it should be time-varying for an orbit)

    R = eye(3);

    d0 = 10; % km

    % System Parameters (Nondim by 1 AU)
    gamma = 10;
    rho = 1E-2;
    [~,Xaug] = ode45(@(t,Xaug) augEOM(t,Xaug,X0,R,Q,d0,gamma,rho,mu,debris_sol),[0,Tf], Xug0, options_ode);
    lambda_f = Xaug(end,7:12);

    % Calc Boundary Condition
    % Costates have to be 0

    % Psi = 0
    Psi = lambda_f;
end