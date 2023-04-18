clc
clear


R_1 = @(a) [1 0 0
            0 cosd(a) sind(a)
            0 -sind(a) cosd(a)];
R_3 = @(a) [cosd(a) sind(a) 0
            -sind(a) cosd(a) 0
            0 0 1];

mu_earth = 3.9860e5;
r_col = [0; 0; 6371*1.2];

a1 = norm(r_col);
a2 = norm(r_col);

v1 = sqrt(mu_earth*(2/norm(r_col) - 1/a1));
v2 = sqrt(mu_earth*(2/norm(r_col) - 1/a2));
v1 = [v1; 0; 0];
% v2 = R_3(-90)*R_1(-90)*R_3(0)*[v2; 0; 0];
v2 = [0;v2;0];

x10 = [r_col; v1];
x20 = [r_col; v2];

opt = odeset("RelTol",1e-8,"AbsTol",1e-8);
[~, x1] = ode45(@(t,x) TwoBodyEOM(t,x,mu_earth), linspace(0,-10000,70000), x10,opt);
[~, x2] = ode45(@(t,x) TwoBodyEOM(t,x,mu_earth), linspace(0,-10000,70000), x20,opt);

r1 = x1(:,1:3);
r2 = x2(:,1:3);

ic1 = x1(end,:);
ic2 = x2(end,:);

figure(1)
plot3(r1(:,1),r1(:,2),r1(:,3))
hold on
plot3(r2(:,1),r2(:,2),r2(:,3))
hold off
axis equal
grid