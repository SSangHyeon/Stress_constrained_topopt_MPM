function [Kt,kp] = G2P(uvw,mpData,nD,quadRho,kp)
penal = 3;
nmp   = length(mpData);                                                     % number of material points
npCnt = 0;                                                                  % counter for the number of entries in Kt
tnSMe = sum([mpData.nSMe]);                                                 % total number of stiffness matrix entries
krow  = zeros(tnSMe,1); kcol=krow; kval=krow;                               % zero the stiffness information

for mp=1:nmp                                                                % material point loop
    nIN = mpData(mp).nIN;                                                   % nodes associated with the material point 
    dNx = mpData(mp).dSvp;                                                  % basis function derivatives (start of lstp)
    nn  = size(dNx,2);                                                      % no. dimensions & no. nodes
    ed  = repmat((nIN-1)*nD,nD,1)+repmat((1:nD).',1,nn);                    % degrees of freedom of nodes (matrix form)
    ed  = reshape(ed,1,nn*nD);                                              % degrees of freedom of nodes (vector form)
    kp{mp,1} = kp{mp,1} .* quadRho(mp)^penal;

    npDoF=(size(ed,1)*size(ed,2))^2;                                        % no. entries in kp
    nnDoF=size(ed,1)*size(ed,2);                                            % no. DoF in kp                        
    krow(npCnt+1:npCnt+npDoF)=repmat(ed',nnDoF,1);                          % row position storage
    kcol(npCnt+1:npCnt+npDoF)=repmat(ed ,nnDoF,1);                          % column position storage
    kval(npCnt+1:npCnt+npDoF)=kp{mp,1};                                           % stiffness storage
    npCnt=npCnt+npDoF;                                                      % number of entries in Kt
end

nDoF=length(uvw);                                                           % number of degrees of freedom
Kt=sparse(krow,kcol,kval,nDoF,nDoF);                                        % form the global stiffness matrix
Kt = (Kt+Kt')/2;
end