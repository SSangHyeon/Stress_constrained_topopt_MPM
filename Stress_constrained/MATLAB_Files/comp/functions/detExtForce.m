function [fext,mpData] = detExtForce(nodes,nD,g,mpData,lstp)

%Global external force determination  
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   23/01/2019
% Description:
% Function to determine the external forces at nodes based on body forces
% and point forces at material points.
%
%--------------------------------------------------------------------------
% [fbdy,mpData] = DETEXTFORCE(coord,etpl,g,eIN,mpData)
%--------------------------------------------------------------------------
% Input(s):
% nodes  - number of nodes (total in mesh)
% nD     - number of dimensions
% g      - gravity
% mpData - material point structured array. Function requires:
%           mpM   : material point mass
%           nIN   : nodes linked to the material point
%           Svp   : basis functions for the material point
%           fp    : point forces at material points
%--------------------------------------------------------------------------
% Ouput(s);
% fext   - external force vector (nodes*nD,1)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

nmp  = size(mpData,2);                                                      % number of material points & dimensions 
fext = zeros(nodes*nD,1);                                                   % zero the external force vector
grav = zeros(nD,1); grav(nD) = -g;                                          % gavity vector
for mp = 1:nmp
   nIN = mpData(mp).nIN;                                                    % nodes associated with MP
   nn  = length(nIN);                                                       % number of nodes influencing the MP
   Svp = mpData(mp).Svp;                                                    % basis functions
%    mpData(mp).fp(1,1) = cos(mpData(mp).ang) * mpData(mp).fp(1,1) - sin(mpData(mp).ang) * mpData(mp).fp(2,1);
%    mpData(mp).fp(2,1) = sin(mpData(mp).ang) * mpData(mp).fp(1,1) + cos(mpData(mp).ang) * mpData(mp).fp(2,1);
%    if lstp >= 2
%         mpData(mp).fp = mpData(mp).R(1:2,1:2) * mpData(mp).fp;
%    end
%    mpData(mp).fp = mpData(mp).F(1:2,1:2) * mpData(mp).fp;
    fp  = (mpData(mp).mpM*grav + mpData(mp).fp)*Svp;                         % material point body & point nodal forces
    ed  = repmat((nIN-1)*nD,nD,1)+repmat((1:nD).',1,nn);                     % nodel degrees of freedom 
   fext(ed) = fext(ed) + fp;                                                % combine into external force vector
end
end