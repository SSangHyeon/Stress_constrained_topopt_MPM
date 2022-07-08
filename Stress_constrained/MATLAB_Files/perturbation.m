clear; clc;
nelx = 160; nely = 80; nelz = 1;
x = 0.5*ones(nely,nelx,nelz);
x = x(:);
d = 1e-4;
e = zeros(nelx*nely,1);
fd = zeros(nelx*nely,1);
pl = 3;
q = 0.5;
p = 8;

for i=1:nelx*nely
    x0 = x;
    x0(i) = x0(i) + d;
    [pnorm,pnorm_sen]=Stress_Sensitivity_Comp(x,pl,q,p);
    [pnorm_f,pnorm_sen]=Stress_Sensitivity_Comp(x0,pl,q,p);
    fd(i) = (pnorm_f-pnorm)/d;
    e(i) = abs((fd(i)-pnorm_sen(i))/pnorm_sen(i));
end