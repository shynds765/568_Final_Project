clear

options_fsolve = optimset('TolFun',1e-8,'TolX',1e-8,'Display','on');
options_ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);

%% Dimensional Parameters
X_debris0 = [1506.88668453313    -11970.5727709491    0    5.37695239639678    6.09571043302428    0].';
X_sc0 = [9146.88869833604    -6365.17383816329    -682.513731171310    0.794523910391561    3.44963774157689    0.369891724256882].';
Tf = 1390.89418662779;

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
R = diag(eye(3));
d0 = 100/l_char; % km

% Collision Cost
gamma = linspace(0.1,2,7);
rho = [1,1E-2,1E-3];

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

figure(1)
plot3(0,0,0)
hold on
plot3(r_debris(:,1),r_debris(:,2),r_debris(:,3),'r-')
hold on
plot3(r_ref(:,1),r_ref(:,2),r_ref(:,3),'k--')
hold on
plot3(r_sc(:,1),r_sc(:,2),r_sc(:,3),'b-')
axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
semilogy(t,dist,'k-')
hold on
yline(d0,'r--')

figure(3)
plot3(0,0,0,'kx')
hold on
plot3(r_sc2debris(:,1),r_sc2debris(:,2),r_sc2debris(:,3),'r-')
hold on
surf(X_ex,Y_ex,Z_ex,'FaceAlpha',0.3, 'FaceColor', 'b')
grid on
axis equal
xlim(2*[-d0,d0])
ylim(2*[-d0,d0])
zlim(2*[-d0,d0])

figure(4)
penalty = gamma(end)/2*(1-tanh((dist-d0)/rho(end)));
plot(t,penalty)
xlabel('Time')
ylabel('Penalty')