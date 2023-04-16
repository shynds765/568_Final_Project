clear

x = linspace(-10,10);

p1 = x.^2+6*x;
p2 = 3*x.^2-8*x;

p_total = p1+p2;

[dom1,ind1] = min(p1);
[dom2,ind2] = min(p2);

[dom_total,ind_total] = min(p_total);

