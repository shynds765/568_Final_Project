function [r0_sat, v0_sat, r0_deb, v0_deb, crashT] = getOrbitalInitialConditions(scenario, wDebris, iDebris, RAANDebris, aDebris, wSat, iSat, RAANSat, aSat)
% Example execution for sample cases 1-5: 
% [r0_sat, v0_sat, r0_deb, v0_deb, crashT] = getOrbitalInitialConditions(1)

% Example execution for user defined case:
% [r0_sat, v0_sat, r0_deb, v0_deb, crashT] = getOrbitalInitialConditions(0, deg2rad(125), deg2rad(92), deg2rad(150), 6786230*100, deg2rad(30), deg2rad(40), deg2rad(15), 6786230*4)

if ~exist('wSat','var')
    wSat = 125*pi/180;
end
if ~exist('iSat','var')
    iSat = 92*pi/180;
end
if ~exist('RAANSat','var')
    RAANSat = 150*pi/180;
end
if ~exist('aSat','var')
    aSat = 6786230*100;
end
if ~exist('w2Debris','var')
    wDebris = 30*pi/180;
end
if ~exist('iDebris','var')
    iDebris = 40*pi/180;
end
if ~exist('RAANDebris','var')
    RAANDebris = 15*pi/180;
end
if ~exist('aDebris','var')
    aDebris = 6786230*4;
end

mu = 3.986*10^5; %GM gravitational constant for Earth
r = 6786230/1000; % [km] semi-major axis ~ 400km orbit 
rinit = [r; 0; 0];
period = 2*pi*sqrt(r^3/mu);
tf = -period/4;
t0 = 0;
t_span = [t0 tf];
v0 = sqrt(mu/r);
vinit = [0; v0; 0];

% Rotation Matrix about the 1st axis
R_1 = @(a) [1 0 0
            0 cos(a) sin(a)
            0 -sin(a) cos(a)];
% Rotation Matrix about the 3rd axis
R_3 = @(a) [cos(a) sin(a) 0
            -sin(a) cos(a) 0
            0 0 1];


switch scenario
    case 0 % user defined inputs
        vElip0 = sqrt(2*mu/r - mu/aDebris);
        vElipInit = [0; vElip0; 0];
        vComet0 = sqrt(2*mu/r - mu/aSat);
        vCometinit = [0; vComet0; 0];
        vCometinit = (R_3(wSat)*R_1(iSat)*R_3(RAANSat)).' * vCometinit;
        vElipInit = (R_3(wDebris)*R_1(iDebris)*R_3(RAANDebris)).' * vElipInit;
        x0 = [rinit(1) rinit(2) rinit(3) vElipInit(1) vElipInit(2) vElipInit(3)];
        xDebris0 = [rinit(1) rinit(2) rinit(3) vCometinit(1) vCometinit(2) vCometinit(3)];
        options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
        [trans,var] = ode45(@TwoBodyEOM, t_span, x0, options, mu);
        [trans2,var2] = ode45(@TwoBodyEOM, t_span, xDebris0, options, mu);

    case 1 % equatorial satellite and equatorial obstacle in opposite directions
        x0 = [rinit(1) rinit(2) rinit(3) vinit(1) vinit(2) vinit(3)];
        xDebris0 = [rinit(1) rinit(2) rinit(3) vinit(1) -vinit(2) vinit(3)];
        options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
        [trans,var] = ode45(@TwoBodyEOM, t_span, x0, options, mu);
        [trans2,var2] = ode45(@TwoBodyEOM, t_span, xDebris0, options, mu);

    case 2 % equatoral satellite and polar obstacle
        x0 = [rinit(1) rinit(2) rinit(3) vinit(1) vinit(2) vinit(3)];
        xDebris0 = [rinit(1) rinit(2) rinit(3) vinit(1) vinit(3) vinit(2)];
        options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
        [trans,var] = ode45(@TwoBodyEOM, t_span, x0, options, mu);
        [trans2,var2] = ode45(@TwoBodyEOM, t_span, xDebris0, options, mu);

    case 3 % equatorial satellite, inclined obstacle
        vDebrisInit = (R_3(wSat)*R_1(iSat)*R_3(RAANSat)).' * vinit;
        x0 = [rinit(1) rinit(2) rinit(3) vinit(1) vinit(2) vinit(3)];
        xDebris0 = [rinit(1) rinit(2) rinit(3) vDebrisInit(1) vDebrisInit(2) vDebrisInit(3)];
        options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
        [trans,var] = ode45(@TwoBodyEOM, t_span, x0, options, mu);
        [trans2,var2] = ode45(@TwoBodyEOM, t_span, xDebris0, options, mu);

    case 4 % equatorial satellite, highly elliptical obstacle
        vComet0 = sqrt(2*mu/r - mu/aSat);
        vCometinit = [0; vComet0; 0];
        vCometinit = (R_3(wSat)*R_1(iSat)*R_3(RAANSat)).' * vCometinit;
        x0 = [rinit(1) rinit(2) rinit(3) vinit(1) vinit(2) vinit(3)];
        xDebris0 = [rinit(1) rinit(2) rinit(3) vCometinit(1) vCometinit(2) vCometinit(3)];
        options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
        [trans,var] = ode45(@TwoBodyEOM, t_span, x0, options, mu);
        [trans2,var2] = ode45(@TwoBodyEOM, t_span, xDebris0, options, mu);

    case 5 % both elliptical
        vElip0 = sqrt(2*mu/r - mu/aDebris);
        vElipInit = [0; vElip0; 0];
        vComet0 = sqrt(2*mu/r - mu/aSat);
        vCometinit = [0; vComet0; 0];
        vCometinit = (R_3(wSat)*R_1(iSat)*R_3(RAANSat)).' * vCometinit;
        x0 = [rinit(1) rinit(2) rinit(3) vElipInit(1) vElipInit(2) vElipInit(3)];
        xDebris0 = [rinit(1) rinit(2) rinit(3) vCometinit(1) vCometinit(2) vCometinit(3)];
        options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
        [trans,var] = ode45(@TwoBodyEOM, t_span, x0, options, mu);
        [trans2,var2] = ode45(@TwoBodyEOM, t_span, xDebris0, options, mu);
end

r0_sat = var(length(var),1:3);
v0_sat = var(length(var),4:6);
r0_deb = var2(length(var2),1:3);
v0_deb = var2(length(var2),4:6);
crashT = -tf;

figure(1)
hold on
plot3(var(:,1), var(:,2), var(:,3))
plot3(var2(:,1), var2(:,2), var2(:,3))
plot3(var(length(var),1), var(length(var),2), var(length(var),3), 'bo')
plot3(var2(length(var2),1), var2(length(var2),2), var2(length(var2),3), 'ro')
axis equal
view(45,45)
legend("Satellite", "Debris")

end