clear

options_fsolve = optimset('TolFun',1e-8,'TolX',1e-8,'Display','on');
options_ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);

mu_earth = 3.9860e5;
mu = 1;
l_char = 6.7862E3; % km
t_char = sqrt(l_char^3/mu_earth);

lambda0_guess = zeros(6,1);
Tf = 2*pi/2;

%% Control Variables
Q = 1*diag(eye(6));

R = diag(eye(3));

d0 = 10/l_char; % km

% Collision Cost
gamma = 1;
rho = 1;

%% Calc Asteroid and Reference Solution
v_c = sqrt(mu_earth/l_char);
v_c_nondim = v_c / (l_char/t_char);

X_debris0 = [-1;0;0;0;v_c_nondim;0];
debris_sol = ode45(@(t,x) TwoBodyEOM(t,x,mu), [0,Tf], X_debris0, options_ode);
t = debris_sol.x;
X_debris = deval(debris_sol,t).';

X0 = [1;0;0;0;v_c_nondim;0];
ref_sol = ode45(@(t,x) TwoBodyEOM(t,x,mu), [0,Tf], X0, options_ode);
X_ref = deval(ref_sol,t).';

%% Run Fsolve
lambda0 = fsolve(@(lamda0) boundaryValueFun(lamda0,Tf,X0,mu,debris_sol,ref_sol,Q,R,d0,gamma,rho),lambda0_guess,options_fsolve);


[~,X] = ode45(@(t,X) augEOM(t,X,Q,R,debris_sol,ref_sol,mu,d0,gamma,rho), t, [X0;lambda0], options_ode);


figure(1)
plot3(X_debris(:,1),X_debris(:,2),X_debris(:,3),'r-')
hold on
plot3(X_ref(:,1),X_ref(:,2),X_ref(:,3),'k--')
hold on
plot3(X(:,1),X(:,2),X(:,3),'b-')
axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')