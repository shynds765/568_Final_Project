function x_dot = TwoBodyEOM(t,x,MU)
r = x(1:3);
v = x(4:6);

% Calculate translational velocity and acceleration
r_dot = v;
v_dot = -MU * r / norm(r)^3;

% Combine State Vector
x_dot = [r_dot; v_dot];
end