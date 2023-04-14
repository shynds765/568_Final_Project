function Xaugdot = augEOM(t,Xaug)
    X = Xaug(1:6);
    lambda = Xaug(7:12);

    x_dot_1body(X(6),nu,u1,u2,u3,x1,x2,x3,x4,x5);

end