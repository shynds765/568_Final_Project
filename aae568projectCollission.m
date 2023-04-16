% Brady Beck 
% AAE 568 
% Project Collission Propogation 
clear all
close all

% Propogating orbit of satellite and rock
a = 6786230; % [m] semi-major axis ~ 400km orbit 
ecc = 0; % eccentricity (deviation of curve from circular)
incl_sat1 = 0; % [deg] inclination of orbit (alternate between 0,180)
incl_sat2 = 0; % [deg] inclination of orbit (alternate between 0,180)
RAAN = 95; % [deg] angle in equatorial plane from xaxis to location of ascending node, not used for equator? 
argp = 93;  % [deg] angle between satellite ascending node and periapsis , not used for circular or equatorial
nu = 0; % [deg] angle between periapsis and current position 
% lon1 = 0:360; % [deg] angle between periapsis and current position 
% lon2 = lon1;
% lon2(1:180) = 180:-1:1;
% lon2(181:361) = 360:-1:180;
lon1 = 0:0.5:180; % [deg] angle between periapsis and current position 
lon2 = 180:-0.5:0;
r_sat1 = zeros(3,181);
r_sat2 = r_sat1; 
for i = 1:size(lon1,2)
    [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl_sat1, RAAN, argp, nu, 'truelon',lon1(i)); % returns 3x1,3x1
    r_sat1(:,i) = r_ijk;
    [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl_sat2, RAAN, argp, nu, 'truelon',lon2(i)); % returns 3x1,3x1
    r_sat2(:,i) = r_ijk;
end

figure(1)
hold on 
grid on
plot3(r_sat1(1,:),r_sat1(2,:),r_sat1(3,:), 'r');
plot3(r_sat1(1,1),r_sat1(2,1),r_sat1(3,1), 'ro', 'MarkerSize',12)
plot3(r_sat2(1,:),r_sat2(2,:),r_sat2(3,:), 'k--');
plot3(r_sat2(1,1),r_sat2(2,1),r_sat2(3,1), 'ko', 'MarkerSize',12)
title('Orbit of Sat1 | Brady Beck')
xlabel('x-direction [m]')
ylabel('y-direction [m]')
zlabel('z-direction [m]')
legend('Orbit1','Start1','Orbit2','Start2')
view(45,45)
hold off

% Plotting relative vector between spacecraft and rock wtih exclusion zone

% Setting up exclusion zone
% Propogating orbit of satellite and rock
exc = 1000; % [m] bound for exclusion zone relative to spacecraft 
rExc = r_sat1 + exc;
sat22sat1 = r_sat2-r_sat1;
exc2sat1 = rExc - r_sat1;
sat221 = zeros(size(lon1));
exc21 = zeros(size(lon1));
exc21_10k = 10000*ones(size(lon1));
for i = 1:size(exc21,2)
    sat221(i) = norm(sat22sat1(:,i));
    exc21(i) = norm(exc2sat1(:,i));
end 

figure(2)
semilogy(lon1,sat221, 'r');
grid on 
hold on 
semilogy(lon1,exc21, 'ko');
semilogy(lon1,exc21_10k, 'bo');
title('Exclusion Zone and Distance to Obstacle, Relative to Sat1 | Brady Beck')
xlabel('index')
ylabel('distance [m]')
legend('Satellite2 to Satellite 1','1km Exclusion to Sattellite', '10k Exclusion to Satellite')
hold off

% Repeating Part 2 for higher order analysis of satellite approaching obstacle
lon1 = 178:0.001:182; % [deg] angle between periapsis and current position 
lon2 = 182:-0.001:178;
r_sat1 = zeros(3,size(lon1,2));
r_sat2 = r_sat1; 
for i = 1:size(lon1,2)
    [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl_sat1, RAAN, argp, nu, 'truelon',lon1(i)); % returns 3x1,3x1
    r_sat1(:,i) = r_ijk;
    [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl_sat2, RAAN, argp, nu, 'truelon',lon2(i)); % returns 3x1,3x1
    r_sat2(:,i) = r_ijk;
end

% Setting up exclusion zone
% Propogating orbit of satellite and rock
exc = 1000; % [m] bound for exclusion zone relative to spacecraft 
rExc = r_sat1 + exc;
sat22sat1 = r_sat2-r_sat1;
[X,Y,Z] = sphere;


% Solving for norm to show 2D plot 
sat221 = zeros(size(lon1));
exc21 = 1000*ones(size(lon1));
exc21_10k = 10000*ones(size(lon1));
for i = 1:size(exc21,2)
    sat221(i) = norm(sat22sat1(:,i));
end 

figure(3)
semilogy(lon1,sat221, 'r');
grid on 
hold on 
semilogy(lon1,exc21, 'ko');
semilogy(lon1,exc21_10k, 'bo');
title('Distance Relative to Sat1, Close View (O: 0.001 degrees) | Brady Beck')
xlabel('index')
ylabel('distance [m]')
legend('Satellite2 to Satellite 1','Exclusion to Sattellite')
hold off

figure(4)
hold on
grid on 
plot3(sat22sat1(1,:),sat22sat1(2,:),sat22sat1(3,:), 'r');
surf(X*10000,Y*10000,Z*10000,'FaceAlpha',0.1, 'FaceColor', 'k')
surf(X*1000,Y*1000,Z*1000,'FaceAlpha',0.3, 'FaceColor', 'b')
title('3D Orbit of Obstacle Relative to Sat1 | Brady Beck')
xlabel('x-direction [m]')
ylabel('y-direction [m]')
zlabel('z-direction [m]')
legend('Satellite 2 to Satellite 1','10 km Exclusion', '1 km Exclusion')
view(45,15)
axis equal
ylim([-20000 20000])
hold off
