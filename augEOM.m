function Xaugdot = augEOM(t,Xaug,Q,R,debris_sol,ref_sol,mu,d0,gamma,rho)
    x = Xaug(1);
    y = Xaug(2);
    z = Xaug(3);
    xdot = Xaug(4);
    ydot = Xaug(5);
    zdot = Xaug(6);
    lambda = Xaug(7:12);

    % Find optimal control
    u = u_star(R(1),R(2),R(3),lambda(4),lambda(5),lambda(6));
    
    % Calculate Position of Debris (in cartesian)
    X_ast = deval(debris_sol,t);
    r_debris = X_ast(1:3);

    % Calculate Reference State
    X_ref = deval(ref_sol,t);
    x_ref = X_ref(1);
    y_ref = X_ref(2);
    z_ref = X_ref(3);
    xdot_ref = X_ref(4);
    ydot_ref = X_ref(5);
    zdot_ref = X_ref(6);

    % State and costate dynamics
    x_dot = EOMwControl(mu,u(1),u(2),u(3),x,xdot,y,ydot,z,zdot);
    lambda_dot = lambdaDot(Q(1),Q(2),Q(3),Q(4),Q(5),Q(6),d0,gamma,lambda(1),lambda(2),lambda(3),lambda(4),lambda(5),lambda(6),...
        mu,r_debris(1),r_debris(2),r_debris(3),rho,x,x_ref,xdot,xdot_ref,y,y_ref,ydot,ydot_ref,z,z_ref,zdot,zdot_ref);

    Xaugdot = [x_dot; lambda_dot];
end