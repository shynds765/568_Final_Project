function [a,e,i] = cartesian2orbitalParams(r,v,mu)
    v_mag = norm(v);
    r_mag = norm(r);

    energy = v_mag^2/2 - mu/r_mag;

    a = -mu/(2*energy);

    h = cross(r,v);
    h_mag = norm(h);

    e = sqrt(1 - h_mag^2/(mu*a));

    h_hat = h/h_mag;

    i = acos(h_hat(3));
end 