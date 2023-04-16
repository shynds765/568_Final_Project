clear

syms mu x y z xdot ydot zdot

r = [x;y;z];
v = [xdot;ydot;zdot];
a = -mu*r/norm(r)^3;

xdot = [v;a];