clear
close all

options_fsolve = optimset('TolFun',1e-8,'TolX',1e-8,'Display','on');
options_ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);

%% Dimensional Parameters
X_debris0 = [1506.88668453313    -11970.5727709491    0    5.37695239639678    6.09571043302428    0].';
X_sc0 = [8510.26572496165    -4584.34548703186    -4990.17150276150    1.29804726770103    2.42612075903804    2.64089142239763].';
% X_debris0 = [-7.2616; 0; -2.3914; -0.0023; 0; 0.0069]*1e3;
% X_sc0 = [0; -7.2616; -2.3914; 0; -0.0023; 0.0069]*1e3;
Tf = 1390.89418662779 * 1.2;

% Characteristic Values
mu_earth = 3.9860e5;
mu = 1;
l_char = norm(X_debris0(1:3)); % km
t_char = sqrt(l_char^3/mu_earth);

%% Nondimensionalize
X_debris0(1:3) = X_debris0(1:3)/l_char;
X_debris0(4:6) = X_debris0(4:6)/(l_char/t_char);

X_sc0(1:3) = X_sc0(1:3)/l_char;
X_sc0(4:6) = X_sc0(4:6)/(l_char/t_char);
Tf = Tf/t_char;

%% Control Variables
Q = 1*diag(eye(6));
R = 1*diag(eye(3));
d0 = 25/l_char; % km

% Collision Cost
% gamma = linspace(0.1,2,7);
% rho = [1,1E-2,1E-3];
gamma = 10;
rho = 1;

%% Calc Asteroid and Reference Solution
debris_sol = ode45(@(t,x) TwoBodyEOM(t,x,mu), [0,Tf], X_debris0, options_ode);
ref_sol = ode45(@(t,x) TwoBodyEOM(t,x,mu), [0,Tf], X_sc0, options_ode);

%% Run Fsolve
lambda0_guess = zeros(6,1);
for i = 1:length(gamma)
    lambda0_guess = fsolve(@(lamda0) boundaryValueFun(lamda0,Tf,X_sc0,mu,debris_sol,ref_sol,Q,R,d0,gamma(i),rho(1)),lambda0_guess,options_fsolve);
end
for i = 1:length(rho)
    lambda0_guess = fsolve(@(lamda0) boundaryValueFun(lamda0,Tf,X_sc0,mu,debris_sol,ref_sol,Q,R,d0,gamma(end),rho(i)),lambda0_guess,options_fsolve);
end
lambda0 = lambda0_guess;

[t,X_sc] = ode45(@(t,X) augEOM(t,X,Q,R,debris_sol,ref_sol,mu,d0,gamma(end),rho(end)), [0,Tf], [X_sc0;lambda0], options_ode);
r_sc = X_sc(:,1:3);

% Find optimal control history
u = zeros(3,length(X_sc));
for i = 1:length(u)
    u(:,i) = u_star(R(1),R(2),R(3),X_sc(i,10),X_sc(i,11),X_sc(i,12))...
        * (l_char/t_char^2);
end

X_debris = deval(debris_sol,t).';
r_debris = X_debris(:,1:3);

X_ref = deval(ref_sol,t).';
r_ref = X_ref(:,1:3);

r_sc2debris = r_debris-r_sc;
dist = vecnorm(r_sc2debris,2,2);

[X_ex,Y_ex,Z_ex] = sphere();
X_ex = d0*X_ex;
Y_ex = d0*Y_ex;
Z_ex = d0*Z_ex;

% Distance from Reference Trajectory
error = zeros(1,length(r_sc));
for i = 1:length(error)
    error(i) = norm(r_sc(i,:) - r_ref(i,:));
end

% Calculating Total Cost from Cost Function Equation
penalty = gamma(end)/2*(1-tanh((dist-d0)/rho(end)));
J = costSum(X_ref-X_sc(:,1:6),Q,u,R,penalty,t);
fprintf('\n-=- Total Cost via Augmented Cost Function -=-\n');
fprintf('Cost with Penalty: %d \n', J(1));
fprintf('Cost without Penalty: %d \n\n', J(2));


figure(1)
plot3(0,0,0)
hold on
plot3(r_debris(:,1)*l_char,r_debris(:,2)*l_char,r_debris(:,3)*l_char,'r-')
hold on
plot3(r_ref(:,1)*l_char,r_ref(:,2)*l_char,r_ref(:,3)*l_char,'k--')
hold on
plot3(r_sc(:,1)*l_char,r_sc(:,2)*l_char,r_sc(:,3)*l_char,'b-')
axis equal
grid on
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

figure(2)
semilogy(t*t_char,dist*l_char,'k-')
hold on
grid on 
xlabel('time [s]')
ylabel('Distance from Satellite [km]')
yline(d0,'r--')

figure(3)
plot3(0,0,0,'kx')
hold on
plot3(r_sc2debris(:,1)*l_char,r_sc2debris(:,2)*l_char,r_sc2debris(:,3)*l_char,'r-')
hold on
surf(X_ex*l_char,Y_ex*l_char,Z_ex*l_char,'FaceAlpha',0.3, 'FaceColor', 'b')
grid on
axis equal
xlim(2*[-d0,d0]*l_char)
ylim(2*[-d0,d0]*l_char)
zlim(2*[-d0,d0]*l_char)
title('Debris Relative to Spacecraft')
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

figure(4)
plot(t,penalty)
xlabel('Time')
ylabel('Penalty')

% Optimal Control History
figure(5)
subplot(3,1,1)
plot(t*t_char,u(1,:))
title("Optimal Control History")
ylabel("x Acceleration [km/s^2]")
xlabel("Time [s]")
grid

subplot(3,1,2)
plot(t*t_char,u(2,:))
ylabel("y Acceleration [km/s^2]")
xlabel("Time [s]")
grid

subplot(3,1,3)
plot(t*t_char,u(3,:))
ylabel("z Acceleration [km/s^2]")
xlabel("Time [s]")
grid

figure(6)
plot(t*t_char,error*l_char)
title("Deviation from Reference Trajectory")
ylabel("Distance from Reference Trajectory [km]")
xlabel("Time [s]")
grid