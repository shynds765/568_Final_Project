clear

options_fsolve = optimset('TolFun',1e-8,'TolX',1e-8,'Display','on');
options_ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);

mu = 3.9860e5;
lambda0_guess = zeros(6,1);
Tf = 2.7818e+03/2;

a0 = 6.7862E3; % km

% In cartesian coords
v_c = sqrt(mu/a0);
X_debris0 = [-a0;0;0;0;v_c;0];
debris_sol = ode45(@(t,x) TwoBodyEOM(t,x,mu), [0,Tf], X_debris0, options_ode);
t = debris_sol.x;
X_debris = deval(debris_sol,t).';

X0 = [a0;0;0;0;v_c;0];
ref_sol = ode45(@(t,x) TwoBodyEOM(t,x,mu), [0,Tf], X0, options_ode);
X_ref = deval(ref_sol,t).';

%lambda0 = fsolve(@(lamda0) boundaryValueFun(lamda0,Tf,X0,mu,debris_sol,ref_sol),lambda0_guess,options_fsolve);

figure(1)
plot3(X_debris(:,1),X_debris(:,2),X_debris(:,3),'r-')
hold on
plot3(X_ref(:,1),X_ref(:,2),X_ref(:,3),'k--')
axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')