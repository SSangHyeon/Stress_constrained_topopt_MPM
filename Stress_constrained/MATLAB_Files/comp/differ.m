clc; clear; close all;
addpath('setup','functions','L_bracket','hole','beam');
[lstps,g,mpData,mesh,mp] = setupGrid_beam_40_20;
penal = 3; volfrac = 0.5; rmin = 1; tol = 1e-9;                                                   % Newton Raphson parameters
[nodes,nD] = size(mesh.coord);                                              % number of nodes and dimensions
[nels,nen] = size(mesh.etpl);                                               % number of elements and nodes/element
nDoF = nodes*nD;                                                            % total number of degrees of freedom
nmp  = length(mpData);                                                      % number of material points  
[mesh,mpData,Refernode,referelem] = elemMPinfo(mesh,mpData);

referelem = referelem';
for i=1:length(referelem)
    referelem(i,2) = mesh.etpl(referelem(i,1),1);
    referelem(i,3) = mesh.etpl(referelem(i,1),2);
    referelem(i,4) = mesh.etpl(referelem(i,1),3);
    referelem(i,5) = mesh.etpl(referelem(i,1),4);
end
iter = 0;
carriers = zeros(3,nmp);
nIN = zeros(nmp,4);

for i = 1:nmp
    carriers(1,i) = mpData(i).mpC(1); carriers(2,i) = mpData(i).mpC(2);
    carriers(3,i) = volfrac;
    nIN(i,:) = mpData(i).nIN; 
end

importfile('ample_beam'); importfile('mesh_beam');  importfile('kp_beam');
importfile('uvw_beam');  importfile('fd_beam'); importfile('fint_beam'); %importfile('kt_5_12800_mp2');
% importfile('BD_5_12800_mp2')
x = volfrac*ones(nmp,1);
K = kp;
  ce = zeros(length(mpData),1);
  change = 1;
  changes = 1e10*ones(10,1);
  edofMat = zeros(nmp,8);
  for i=1:nmp
      J(i) = length(mpData(i).nIN);
      for j=1:J(i)
        edofMat(i,2*j-1) = mpData(i).nIN(j)*2-1; edofMat(i,2*j) = mpData(i).nIN(j)*2;
      end
  end
den = carriers(3,:);

%% Iteration Loop
while change > 1e-3 && iter < 100
  iter = iter + 1; 
  U = zeros(nDoF,1);
  [Kt,kp] = G2P(tot_uvw,mpData,nD,den',K);
  U(fd)=Kt(fd,fd)\fint(fd);
  for i=1:nmp
      mat = edofMat(i,1:2*J(1,i));
      ce(i) = sum(U(mat)'*kp{i,1}*U(mat));
  end
%   c = sum(ce.*penal.*(tol+(1-den').^penal*(1-tol)));
  c = sum(ce);
  v = sum((den)/numel(den));
  volume = sum(carriers(3,:)/numel(carriers(3,:)));
  dCdRhoNaive = -penal*((1-den').^(penal-1).*ce);
  max_sen = max(abs(dCdRhoNaive));
  sensitivity = abs(dCdRhoNaive./max_sen);
  dC_norm = max(abs(dCdRhoNaive));

  %% 99-line code sensitivity analysis & OC check
  [ce] = check(carriers(1:2,:)',den',ce,rmin);
  xnew = OC(nmp,den',volfrac,ce);
  den = xnew';
  %%
  changes(mod(iter,10)+1) = c;
  change = min(1,(max(changes)-min(changes))/min(changes));
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f dC_norm: %7.3f relCh: %7.3f\n',iter,c,v, dC_norm, change);
  %% 
  figure(1)
  scatter(carriers(1,:),carriers(2,:),5,(1-den),"filled",'o');
  colormap('gray');
  caxis([0 1]); axis equal; axis off; drawnow;
  den_total(:,iter) = 1-den;
  obj(iter) = c;
  vol(iter) = volume;
  volume_quad(iter) = v;
  iteration(iter) = iter;

%   if iter == 1
%     figure(2)
%     scatter(carriers(1,:),carriers(2,:),40,(sensitivity),"filled","s"); colormap('turbo'); axis equal; axis off; drawnow
%   end
end
figure(3)
plot(iteration,obj,'b-o');
xlabel('iteration'); ylabel('compliance'); grid on; 
hold on;
% plot(iteration,vol,'r-x');
% xlabel('iteration'); ylabel('Total volume of carrier'); grid on;
% ylim([0 1])
figure(5)
for i=1:nmp
    x = abs(mpData(i).sig(1)); y = abs(mpData(i).sig(2)); z = abs(mpData(i).sig(3)); xy = abs(mpData(i).sig(4));
    von(i) = sqrt((((x-y)^2)+((y-z)^2)+((z-x)^2)+(6*(xy^2)))/2);
    von(i) = sqrt(x^2 - x*z +z^2 +3*y^2);
end
% pnorm = (sum(von.^penal))^(1/p);
scatter(carriers(1,:),carriers(2,:),20,von,'o',"filled");
colormap('jet')
axis equal
axis off
% colorbar('southoutside');
% figure(5)
% scatter(carriers(1,:),carriers(2,:),3,'filled','o','k'), colormap('gray'); caxis([0 1]); axis equal
% axis off
% figure(6)
% scatter(carriers(1,:),carriers(2,:),20,'filled','o'); colormap('gray'); caxis([0 1]); axis equal; axis off; drawnow;
% for i=1:iter
%     scatter(carriers(1,:),carriers(2,:),20,(den_total(:,i)),"filled",'o'); colormap('gray'); caxis([0 1]); axis equal; axis off; drawnow;
%     pause(0.1);
% end