clc
clear

syms gamma d0 nu n rho 
syms pos_obj [3 1]
syms x [6 1]
syms x_target [6 1]
syms u [3 1]
syms Q1 Q2 Q3 Q4 Q5 Q6 
syms R1 R2 R3
syms lambda [6 1]

Q = diag([Q1,Q2,Q3,Q4,Q5,Q6]);
R = diag([R1,R2,R3]);

mu = 3.9860e5;
deltax = x - x_target;

% Rotation Matrix about the 1st axis
R1 = @(a) [1 0 0
            0 cos(a) sin(a)
            0 -sin(a) cos(a)];
% Rotation Matrix about the 3rd axis
R3 = @(a) [cos(a) sin(a) 0
            -sin(a) cos(a) 0
            0 0 1];

% Express distance between spacecraft in terms of states
a = x(1);
e = x(2);
i = x(3);
RAAN = x(4);
w = x(5);
M = x(6);

p = a*(1-e^2);

r = p/(1+e*cos(nu));
pos_sc = [r*cos(nu);
          r*sin(nu);
          0];  % Positision in 2d
pos_sc = (R3(w)*R1(i)*R3(RAAN)).' * pos_sc; % 3d position

r_c = norm(pos_obj-pos_sc); % Distance between spacecrafts

% Cost function
Pc = gamma/2*(1-tanh((r_c-d0)/rho)); % Exclusion Cost
Pp = deltax.' * Q * deltax; % Tracking Cost
L = u.' * R * u + Pc + Pp;

% Dynamics
r = norm(pos_sc);
h = sqrt(mu*p);
b = a*sqrt(1-e^2);
B = (1/h) * [2*a^2*e*sin(nu) 2*a^2*p/r 0
            p*sin(nu) (p+r)*cos(nu)+r*e 0
            0 0 r*cos(nu + w)
            0 0 (r*sin(nu + w))/(sin(i))
            -p*cos(nu)/e (p+r)*sin(nu)/e -r*sin(nu+w)/tan(i)
            (b*p*cos(nu)/(a*e))-(2*b*r/a) -b*(p+r)*sin(nu)/(a*e) 0]; 
x_dot = [zeros(5,1); n] + B*u;

% Hamiltonian
H = L + lambda.' * x_dot;

% Costate Dynamics
lambda_dot = -gradient(H,x1);

% Optimal Control 
u_star = -R^(-1)*B.'*lambda;

%matlabFunction(lambda_dot,'File','lambda_dot_1body');
%matlabFunction(lambda_dot,'File','x_dot_1body'); 
%matlabFunction(u_star,'File','u_star_1body'); 