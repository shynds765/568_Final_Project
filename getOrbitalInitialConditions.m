function [r0_sat, v0_sat, r0_debris, v0_debris] = getOrbitalInitialConditions(scenario, w, i, RAAN)

if ~exist('w','var')
    w = 30*pi/180;
end
if ~exist('i','var')
    i = 5*pi/180;
end
if ~exist('RAAN2','var')
    RAAN = 15*pi/180;
end

mu = 3.986*10^5; %GM gravitational constant for Earth
r = 6786230/1000; % [m] semi-major axis ~ 400km orbit 
rinit = [r; 0; 0];
r_dot0 = 0;
period = 2*pi*sqrt(r^3/mu);
tf = -period/4;
t0 = 0;
t_span = [t0 tf];
v0 = sqrt(mu/r);
vinit = [0; v0; 0];
theta0 = 0;
theta_dot0 = v0 / r;
a = 6786230*10;
vComet0 = sqrt(2*mu/r - mu/a);
vCometinit = [0; vComet0; 0];

% Rotation Matrix about the 1st axis
R_1 = @(a) [1 0 0
            0 cos(a) sin(a)
            0 -sin(a) cos(a)];
% Rotation Matrix about the 3rd axis
R_3 = @(a) [cos(a) sin(a) 0
            -sin(a) cos(a) 0
            0 0 1];


switch scenario
    case 0
    case 1 % equatorial satellite and equatorial obstacle in opposite directions
    case 2 % equatoral satellite and polar obstacle
    case 3 % equatorial satellite, inclined obstacle
    case 4 % equatorial satellite, highly elliptical obstacle
    case 5 % both elliptical



vinit = (R_3(w)*R_1(i)*R_3(RAAN)).' * vinit;
x0 = [rinit(1) rinit(2) rinit(3) vinit(1) vinit(2) vinit(3)];
xComet0 = [rinit(1) rinit(2) rinit(3) vCometinit(1) vCometinit(2) vCometinit(3)];
options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
[trans,var] = ode45(@TwoBodyEOM, t_span, x0, options, mu);
[trans2,var2] = ode45(@TwoBodyEOM, t_span, xComet0, options, mu);

figure(1)
hold on
plot3(var(:,1), var(:,2), var(:,3))
plot3(var2(:,1), var2(:,2), var2(:,3))
plot3(var(length(var),1), var(length(var),2), var(length(var),3), 'bo')
plot3(var2(length(var2),1), var2(length(var2),2), var2(length(var2),3), 'ro')
axis equal
view(45,45)



end