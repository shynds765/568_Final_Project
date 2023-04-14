function Xaugdot = augEOM(t,Xaug, R)
    x = Xaug(1:6);
    lambda = Xaug(7:12);

    u = u_star_1body(R(1),R(2),R(3),lambda(1),lambda(2),lambda(3),lambda(4),lambda(5),lambda(6),x(1),x(2),x(3),x(4),x(5),x(6));
    x_dot = x_dot_1body(x(6),nu,u(1),u(2),u(3),x(1),x(2),x(3),x(4),x(5));

    Xaugdot = [x_dot(1);
               x_dot(2);
               x_dot(3);
               x_dot(4);
               x_dot(5);
               x_dot(6);
               lambda(1);
               lambda(2);
               lambda(3);
               lambda(4);
               lambda(5);
               lambda(6)];
end