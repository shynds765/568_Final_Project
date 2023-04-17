% Brady Beck 
% AAE 568 
% Project Collission Propogation 
clear all
close all

% Propogating orbit of satellite and obstacle
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
lon1 = 0:0.5:360; % [deg] angle between periapsis and current position 
lon2 = 180:-0.5:0;
r_sat1 = zeros(3,181);
r_sat2 = r_sat1; 
v_sat1 = r_sat1;
v_sat2 = r_sat1;
for i = 1:size(lon1,2)
    [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl_sat1, RAAN, argp, nu, 'truelon',lon1(i)); % returns 3x1,3x1
    r_sat1(:,i) = r_ijk;
    v_sat1(:,i) = v_ijk;
%     [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl_sat2, RAAN, argp, nu, 'truelon',lon2(i)); % returns 3x1,3x1
%     r_sat2(:,i) = r_ijk;
%     v_sat2(:,i) = v_ijk;
end

% Rotation Matrix about the 1st axis
R_1 = @(a) [1 0 0
            0 cos(a) sin(a)
            0 -sin(a) cos(a)];
% Rotation Matrix about the 3rd axis
R_3 = @(a) [cos(a) sin(a) 0
            -sin(a) cos(a) 0
            0 0 1];
w = 30*pi/180;
i = 5*pi/180;
RAAN = 15*pi/180;
rInc = (R_3(w)*R_1(i)*R_3(RAAN)).' * r_sat1; 
vInc = (R_3(w)*R_1(i)*R_3(RAAN)).' * v_sat1;

figure(1)
hold on 
grid on
plot3(r_sat1(1,:),r_sat1(2,:),r_sat1(3,:), 'r');
plot3(rInc(1,:),rInc(2,:),rInc(3,:), 'b');
title('Orbit of Sat1 | Brady Beck')
xlabel('x-direction [m]')
ylabel('y-direction [m]')
zlabel('z-direction [m]')
legend('Planar', 'Inclined')
view(45,45)
hold off
fprintf('\n-=- Orbit Conditions [deg] -=-\n');
fprintf('Argument of Perigee: %d\n', rad2deg(w));
fprintf('Inclination: %d\n', rad2deg(i));
fprintf('RAAN: %d\n', rad2deg(RAAN));
fprintf('\n-=- Initial Conditions -=-\n');
fprintf('r: ');
disp(rInc(:,1)');
fprintf('v: ');
disp(vInc(:,1)');