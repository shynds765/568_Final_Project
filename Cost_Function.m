clc
clear

syms gamma d0 nu
syms x1 x2 [6 1] matrix
syms u [12 1] matrix
syms Q1 Q2 [6 6] matrix
syms R [12 12] matrix
syms lamda [12 1]

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
r1 = [(p1 * cos(nu)) / (1 + x1(2)*cos(nu))
    (p1 * sin(nu)) / (1 + x1(2)*cos(nu))
    0];  % Positision in 2d
r1 = R3(-x1(4))*R1(-x1(3))*R3(-x1(5)) * r1; % 3d position
p2 = x2(1) * (1 - x2(2)^2);  % Orbital Parameter
r2 = [(p2 * cos(nu)) / (1 + x2(2)*cos(nu))
    (p2 * sin(nu)) / (1 + x2(2)*cos(nu))
    0];  % Positision in 2d
r2 = R3(-x2(4))*R1(-x2(3))*R3(-x2(5)) * r2; % 3d position
r = norm(r2-r1); % Distance between spacecrafts

% Cost function
Pc = gamma * (1 - (r/d0)^2)^3;
Pp1 = -x1.' * Q1 * x1;
Pp2 = -x2.' * Q2 * x2;
J = u.' * R * u + Pc + Pp1 + Pp2;