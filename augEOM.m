function Xaugdot = augEOM(t,Xaug,R,Q,d0,gamma,rho)
    x = Xaug(1:6);
    lambda = Xaug(7:12);

    % Find optimal control
    u = u_star_1body(R(1,1),R(2,2),R(3,3),lambda(1),lambda(2),lambda(3),...
        lambda(4),lambda(5),lambda(6),x(1),x(2),x(3),x(4),x(5),x(6));
    
    % State and costate dynamics
    x_dot = x_dot_1body(x(6),u(1),u(2),u(3),x(1),x(2),x(3),x(4),x(5),x(6));
    lambda_dot = lambda_dot_1body(Q(1,1),Q(2,2),Q(3,3),Q(4,4),Q(5,5),...
        Q(6,6),d0,gamma,lambda(1),lambda(2),lambda(3),lambda(4),...
        lambda(5),lambda(6),pos_obj1,pos_obj2,pos_obj3,rho,u(1),u(2),...
        u(3),x(1),x(2),x(3),x(4),x(5),x(6),x_target1,x_target2,...
        x_target3,x_target4,x_target5,x_target6);

    Xaugdot = [x_dot; lambda_dot];
end