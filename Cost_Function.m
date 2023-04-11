clc
clear

syms gamma d0 nu
syms x1 x2 [6 1]
syms u [6 1]
syms Q1 Q2 R [6 6]
syms lambda [12 1]

% Rotation Matrix about the 1st axis
R1 = @(a) [1 0 0
    0 cos(a) sin(a)
    0 -sin(a) cos(a)];
% Rotation Matrix about the 3rd axis
R3 = @(a) [cos(a) sin(a) 0
    -sin(a) cos(a) 0
    0 0 1];

% Express distance between spacecraft in terms of states
p1 = x1(1) * (1 - x1(2)^2);  % Orbital Parameter
pos1 = [(p1 * cos(nu)) / (1 + x1(2)*cos(nu))
    (p1 * sin(nu)) / (1 + x1(2)*cos(nu))
    0];  % Positision in 2d
pos1 = R3(-x1(4))*R1(-x1(3))*R3(-x1(5)) * pos1; % 3d position
p2 = x2(1) * (1 - x2(2)^2);  % Orbital Parameter
pos2 = [(p2 * cos(nu)) / (1 + x2(2)*cos(nu))
    (p2 * sin(nu)) / (1 + x2(2)*cos(nu))
    0];  % Positision in 2d
pos2 = R3(-x2(4))*R1(-x2(3))*R3(-x2(5)) * pos2; % 3d position
r = norm(pos2-pos1); % Distance between spacecrafts

% Cost function
Pc = gamma * (1 - (r/d0)^2)^3;
Pp1 = -x1.' * Q1 * x1;
Pp2 = -x2.' * Q2 * x2;
J = u.' * R * u + Pc + Pp1 + Pp2;

% Dynamics
a1 = x1(1);
e1 = x1(2);
r1 = norm(pos1);
omega1 = x1(5);
i1 = x1(3);
b1 = a1 * (1 - e1^2).^(.5);
a2 = x2(1);
e2 = x2(2);
r2 = norm(pos2);
omega2 = x2(5);
i2 = x2(3);
b2 = a2 * (1 - e2^2).^(.5);
B1 = [2*a1^2*e1*sin(nu) 2*a1^2*p1/r1 0
    p1*sin(nu) (p1+r1)*cos(nu)+r1*e1 0
    0 0 r1*cos(nu + omega1)
    0 0 (r1*sin(nu + omega1))/(sin(i1))
    -p1*cos(nu)/e1 (p1+r1)*sin(nu)/e1 -r1*sin(nu+omega1)/tan(i1)
    b1*p1*cos(nu)/(a1*e1)-2*b1*r1/a1 -b1*(p1+r1)*sin(nu)/(a1*e1) 0]; % still need 1/h
B2 = [2*a2^2*e2*sin(nu) 2*a2^2*p2/r2 0
    p2*sin(nu) (p2+r2)*cos(nu)+r2*e2 0
    0 0 r2*cos(nu + omega2)
    0 0 (r2*sin(nu + omega2))/(sin(i2))
    -p2*cos(nu)/e2 (p2+r2)*sin(nu)/e2 -r2*sin(nu+omega2)/tan(i2)
    b2*p2*cos(nu)/(a2*e2)-2*b2*r2/a2 -b2*(p2+r2)*sin(nu)/(a2*e2) 0]; % still need 1/h
B = blkdiag(B1,B2);
x_dot = B*u;

% Hamiltonian
H = J + lambda.' * x_dot;

% Costate Dynamics
lambda_dot = -[gradient(H,x1); gradient(H,x2)];

% Optimal control
dHdu = gradient(H,u);