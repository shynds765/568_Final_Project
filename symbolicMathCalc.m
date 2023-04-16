clear

syms mu x y z xdot ydot zdot
syms x_ref y_ref z_ref xdot_ref ydot_ref zdot_ref
syms u1 u2 u3
syms Q1 Q2 Q3 Q4 Q5 Q6 
syms R1 R2 R3
syms lambda [6 1]
syms gamma rho d0
syms r_debris [3 1]


Q = diag([Q1,Q2,Q3,Q4,Q5,Q6]);
R = diag([R1,R2,R3]);

X_target = [x_ref;y_ref;z_ref;xdot_ref;ydot_ref;zdot_ref]; 

r = [x;y;z];
v = [xdot;ydot;zdot];
X = [r;v];
a = -mu*r/norm(r)^3;

Xdot = [v;a];

u = [u1;u2;u3];
B = [zeros(3);eye(3)];

f = Xdot + B*u;

% Cost function
deltaX = X - X_target;
r_c = norm(r - r_debris); % Distance between spacecrafts
Pc = gamma/2*(1-tanh((r_c-d0)/rho)); % Exclusion Cost
Pp = deltaX.' * Q * deltaX; % Tracking Cost
L = .5 * (u.' * R * u + Pp) + Pc;

% Hamiltonian
H = L + lambda.' * Xdot;

% Costate Dynamics
lambda_dot = -gradient(H,X);

% Optimal Control 
u_star = -R^(-1)*B.'*lambda;

%matlabFunction(lambda_dot,'File','lambdaDot');
%matlabFunction(f,'File','EOMwControl'); 
%matlabFunction(u_star,'File','u_star'); 