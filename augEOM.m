function Xaugdot = augEOM(t,Xaug,x_target,R,Q,d0,gamma,rho,mu,debris_sol)
    x = Xaug(1:6);
    lambda = Xaug(7:12);

    % Find optimal control
    u = u_star_1body(R(1,1),R(2,2),R(3,3),lambda(1),lambda(2),lambda(3),...
        lambda(4),lambda(5),lambda(6),x(1),x(2),x(3),x(4),x(5),x(6));
    
    % Calculate Position of Debris (in cartesian)
    X_ast = deval(debris_sol,t);
    pos_obj = X_ast(1:3);

    n = sqrt(mu/x(1)^3);

    % State and costate dynamics
    x_dot = x_dot_1body(n,u(1),u(2),u(3),x(1),x(2),x(3),x(4),x(5),x(6));
    lambda_dot = lambda_dot_1body(Q(1,1),Q(2,2),Q(3,3),Q(4,4),Q(5,5),...
        Q(6,6),d0,gamma,lambda(1),lambda(2),lambda(3),lambda(4),...
        lambda(5),lambda(6),pos_obj(1),pos_obj(2),pos_obj(3),rho,u(1),u(2),...
        u(3),x(1),x(2),x(3),x(4),x(5),x(6),x_target(1),x_target(2),...
        x_target(3),x_target(4),x_target(5),x_target(6));

    Xaugdot = [x_dot; lambda_dot];
end