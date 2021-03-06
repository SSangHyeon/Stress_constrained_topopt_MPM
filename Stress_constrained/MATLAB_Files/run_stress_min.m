clc; clear; close('all'); warning('off','all')
addpath('setup','functions','comp','MMA');

importfile('mpdata_L');
nmp  = length(mpData);
carriers = zeros(3,nmp);

for i = 1:nmp
    carriers(1,i) = mpData(i).mpC(1); carriers(2,i) = mpData(i).mpC(2);
end

importfile('ample_L'); importfile('mesh_L');  importfile('kp_L');importfile('uvw_L');
importfile('fd_L'); importfile('fint_L'); importfile('kt_L'); importfile('B_L'); importfile('D_L');

rmin =1;
[Hs,H]=prepare_filter(rmin,carriers);

x=0.5*ones(nmp,1);
m =2; %% 1 = volume constraint, x = clustering stress constraint
epsimin = 0.0000001;
n=length(x(:));
xval=x(:);
xold1   = xval;
xold2   = xval;
xlb = 1e-3*ones(n,1);
xub = 1*ones(n,1);
xmin    = xlb;
xmax    = xub;
low     = xlb;
upp     = xub;
c       = 1e3*ones(m,1);
d       = [0]';
a0      = 1;
a       = zeros(m,1);
raa0    = 0.0001;
raa     = 0.0001;
raa0eps = 0.0000001;
raaeps  = 0.0000001;
outeriter = 0;
maxoutit  =100;
kkttol  = 0;
x_his=zeros(nmp,maxoutit);
if outeriter < 0.5
[f0val,df0dx,fval,dfdx,pnorm]=stress_minimize(xval,B,D,carriers,fd,kp,mpData,tot_uvw,Hs,H,fint);
innerit=0;
outvector1 = [outeriter innerit xval'];
outvector2 = [f0val fval'];
end
kktnorm = kkttol+1;
outit = 0;
f = zeros(maxoutit,4);

while  outit < maxoutit 
outit   = outit+1;
outeriter = outeriter+1;
%% GCMMA %%
%%%% The parameters low, upp, raa0 and raa are calculated:
[low,upp,raa0,raa] = ...
asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
raa0,raa,raa0eps,raaeps,df0dx',dfdx);
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
xold2 = xold1;
xold1 = xval;
xval  = xmma;
%% Update %%
[f0val,df0dx,fval,dfdx,pnorm]=stress_minimize(xval,B,D,carriers,fd,kp,mpData,tot_uvw,Hs,H,fint);

%% PRINT RESULTS
fprintf(' It.:%5i      Compliance.:%11.4f   Vol.:%7.3f   P-norm.:%7.3f  \n',outit,f0val,mean(xval(:)),pnorm);
f(outit,1) = outit; f(outit,2) = f0val; f(outit,3) = mean(xval(:)); f(outit,4) = pnorm;
%%%% The residual vector of the KKT conditions is calculated:
[residu,kktnorm,residumax] = ...
kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
outvector1 = [outeriter innerit xval'];
outvector2 = [f0val fval'];
x_his(:,outit)=xmma;
end
%% Wang interpolation
% beta1 = 1; 
% rho0  = 0.01*ones(nmp,1); 
% gamma = tanh(beta1*rho0) + tanh(beta1 * (xval - rho0));
% xval = gamma ./(tanh(beta1*rho0) + tanh(beta1*(1.0-rho0)));
% figure(7); scatter(carriers(1,:),carriers(2,:),40,-xval,"filled","s");axis equal;axis off; colormap('gray');title('material layout');drawnow;
%% figure
figure(6)
yyaxis left
plot(f(1:outit,1),f(1:outit,2),'-o'); ylabel('Compliance & P-norm stress/100');
hold on;
plot(f(1:outit,1),f(1:outit,4)./100,'-*');
yyaxis right
plot(f(1:outit,1),f(1:outit,3),'-x'); ylabel('volume');
xlabel('iteration');
legend('compliance','P-norm stress','volume');