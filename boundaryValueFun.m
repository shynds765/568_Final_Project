function Psi = boundaryValueFun(lamda0,Tf,X0)
    options_ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);
    Xug0 = [X0;lamda0];

    % System Parameters (Nondim by 1 AU)
    [~,Xaug] = ode45(@(t,Xaug) augEOM(t,Xaug,minimize_case),[0,Tf], Xug0, options_ode);

    % Calc Boundary Condition
    % lambda = 0

    % Psi = 0
    Psi = [r_Mars;v_Mars] - state_f;
end