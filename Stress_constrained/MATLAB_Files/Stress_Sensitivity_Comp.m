function [pnorm,pnorm_sen,von,x,ce]=Stress_Sensitivity_Comp(x,penal,q,p,B,D,fd,kp,mpData,tot_uvw,fint)
addpath('setup','functions');
nmp  = length(mpData);                                                      % number of material points  
von = zeros(nmp,1);
edofMat = zeros(nmp,8);
  for i=1:nmp
      J = length(mpData(i).nIN);
      for j=1:J
        edofMat(i,2*j-1) = mpData(i).nIN(j)*2-1; edofMat(i,2*j) = mpData(i).nIN(j)*2;
      end
  end
S = zeros(nmp,3);
U = zeros(length(tot_uvw),1);
ce = zeros(length(mpData),1);
[Kt,kp] = G2P(tot_uvw,mpData,2,x,kp);
U(fd)=Kt(fd,fd)\fint(fd);
for i=1:nmp
    mat = edofMat(i,:);
    ce(i) = U(mat)'*kp{i,1}*U(mat);
end
for i=1:nmp
    temp= (x(i)^q)*(D(:,:,1)*B(:,:,i)*U(edofMat(i,:)))';
    S(i,:) = temp;
    von(i) =sqrt(temp(1)^2 + temp(2)^2 - temp(1)*temp(2) +3*temp(3)^2);
end
DvmDs=zeros(nmp,3);
for i=1:nmp
    DvmDs(i,1)=(1/(2*von(i)))*(2*S(i,1)-S(i,2));
    DvmDs(i,2)=(1/(2*von(i)))*(2*S(i,2)-S(i,1));
    DvmDs(i,3)=(3/von(i))*S(i,3); 
end
dpn_dvms=(sum(von.^p))^(1/p-1);
index_matrix=edofMat';
pnorm=(sum(von.^p))^(1/p);
beta=zeros(nmp,1);
for i=1:nmp
    u=reshape(U(edofMat(i,:),:)',[],1);
    beta(i)=q*(x(i)^(q-1)) * von(i)^(p-1) * DvmDs(i,:)*D(:,:,i)*B(:,:,i)*u;
end
T1=dpn_dvms*beta;
gama=zeros(nmp*2,1);
for i=1:nmp
    index=index_matrix(:,i);
    gama(index)=gama(index)+ x(i)^q * dpn_dvms*B(:,:,i)'*D(:,:,i)'*DvmDs(i,:)'* von(i).^(p-1);
end
lamda=zeros(nmp*2,1);
lamda(fd,:)=Kt(fd,fd)\gama(fd,:);
T2=zeros(nmp,1);
for i=1:nmp
    index=index_matrix(:,i);
    T2(i)=-lamda(index)'*penal* x(i)^(penal-1) *kp{i,1}*U(index);
end
pnorm_sen=T1+T2;
end