function [f0val,df0dx,fval,dfdx,pnorm]=stress_minimize(x,B,D,carriers,Fd,kp,mpData,uvw,Hs,H,fint)
nmp = length(carriers);
pl=3;q=0.5; p=40;
x(:)=(H*x(:))./Hs;
[pnorm,pnorm_sen,MISES,x,ce]=Stress_Sensitivity_Comp(x,pl,q,p,B,D,Fd,kp,mpData,uvw,fint);
% x(:)=H*(x(:)./Hs);
pnorm_sen = H*(pnorm_sen(:)./Hs);
ce = H*(ce(:)./Hs);
% meaning = mean(abs(pnorm_sen));
% cons = find(abs(pnorm_sen)>mean(abs(pnorm_sen)*100));
% for i=1:length(cons)
%     pnorm_sen(cons(i)) = sign(pnorm_sen(cons(i)))*meaning*100;
% end
% for i=1:nmp
%     pnorm_sen(i) = sign(pnorm_sen(i))*log(abs(pnorm_sen(i)));
% end

figure(1); scatter(carriers(1,:),carriers(2,:),40,-x,"filled","s");axis equal;axis off; colormap('gray');title('material layout');drawnow;
figure(2); scatter(carriers(1,:),carriers(2,:),40,MISES,"filled","s"); axis equal;axis off; colormap('jet');title('Von-Mises Stress');colorbar;drawnow;
% figure(3); scatter(carriers(1,:),carriers(2,:),40,pnorm_sen,"filled","s"); axis equal;axis off; colormap('jet');title('Sensitivity');colorbar;drawnow;
% figure(4); plot3(carriers(1,:),carriers(2,:),pnorm_sen,'.'); colormap('jet'); xlabel('x'); ylabel('y'); zlabel('sensitivity');drawnow;
dv = ones(nmp,1)/(nmp);
dv(:) = (H*dv(:))./Hs;
fval=[mean(x)-0.5;pnorm/10];
dfdx=[dv(:)';pnorm_sen'];
df0dx=ce;
f0val=sum(ce);