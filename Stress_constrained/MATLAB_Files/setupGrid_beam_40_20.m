function [lstps,g,mpData,mesh,mp] = setupGrid_beam_40_20
%% Analysis parameters
E=12e6;   v=0.3;   fc=4e3;                                                 % Young's modulus, Poisson's ratio, yield strength   
mCst=[E v fc];                                                              % material constants
g=0;                                                                       % gravity
rho=1;                                                                      % material density
P=-5e3;                                                                   % applied end load
lstps=10;                                                                   % number of loadsteps
a = 1;                                                                      % element multiplier
nelsx=42;                                                                 % number of elements in the x direction
nelsy=40;                                                                 % number of elements in the y direction
ly=40;  lx=42;                                                              % domain dimensions
d=20;  l=40;                                                                 % beam dimensions
mp=4;                                                                       % number of material points in each direction per element
mpType = 1;                                                                 % material point type: 1 = MPM, 2 = GIMP
cmType = 1;                                                                 % constitutive model: 1 = elastic, 2 = vM plasticity
mpSte=0;
%% Mesh generation
[etpl,coord] = formCoord2D(nelsx,nelsy,lx,ly);                              % background mesh generation
[~,nen]      = size(etpl);                                                  % number of element nodes
[nodes,nD]   = size(coord);                                                 % number of nodes and dimensions
h            = [lx ly]./[nelsx nelsy];                                      % element lengths in each direction

%% Boundary conditions on background mesh
bc = zeros(nodes*nD,2);                                                     % generate empty bc matrix
for node=1:nodes                                                            % loop over nodes
  if coord(node,1)==0                                                       % roller (x=0)
    bc(node*2-1,:)=[node*2-1 0];
    bc(node*2  ,:)=[node*2   0];
  end
%   if coord(node,1)==0 && coord(node,2)==(ly-d/2)                            % mid-depth pin
%     bc(node*2  ,:)=[node*2   0];
%   end
end
bc = bc(bc(:,1)>0,:);                                                       % remove empty part of bc

%% Mesh data structure generation
mesh.etpl  = etpl;                                                          % element topology
mesh.coord = coord;                                                         % nodal coordinates
mesh.bc    = bc;                                                            % boundary conditions
mesh.h     = h;                                                             % mesh size

%% Material point generation
ngp    = mp^nD;                                                             % number of material points per element
GpLoc  = detMpPos(mp,nD);                                                   % local MP locations
N      = shapefunc(nen,GpLoc,nD);                                           % basis functions for the material points
[etplmp,coordmp] = formCoord2D(40,20,l,d);                               % mesh for MP generation
coordmp(:,2)=coordmp(:,2)+(ly-d);                                           % adjust MP locations (vertical)
nelsmp = size(etplmp,1);                                                    % no. elements populated with material points
nmp    = ngp*nelsmp;                                                        % total number of mterial points

mpC=zeros(nmp,nD);                                                          % zero MP coordinates
for nel=1:nelsmp
  indx=(nel-1)*ngp+1:nel*ngp;                                               % MP locations within mpC
  eC=coordmp(etplmp(nel,:),:);                                              % element coordinates
  mpPos=N*eC;                                                               % global MP coordinates
  mpC(indx,:)=mpPos;                                                        % store MP positions
end
lp = zeros(nmp,2);                                                          % zero domain lengths
lp(:,1) = h(1)/(2*mp);                                                      % domain half length x-direction
lp(:,2) = h(2)/(2*mp);                                                      % domain half length y-direction
vp      = 2^nD*lp(:,1).*lp(:,2);                                            % volume associated with each material point
count =0;
%% Material point structure generation
for mp = nmp:-1:1                                                           % loop backwards over MPs so array doesn't change size
  mpData(mp).mpType = mpType;                                               % material point type: 1 = MPM, 2 = GIMP
  mpData(mp).cmType = cmType;                                               % constitutive model: 1 = elastic, 2 = vM plasticity
  mpData(mp).mpC    = mpC(mp,:);                                            % material point coordinates
  mpData(mp).vp     = vp(mp);                                               % material point volume
  mpData(mp).vp0    = vp(mp);                                               % material point initial volume
  mpData(mp).mpM    = vp(mp)*rho;                                           % material point mass
  mpData(mp).nIN    = zeros(nen,1);                                         % nodes associated with the material point
  mpData(mp).eIN    = 0;                                                    % element associated with the material point
  mpData(mp).Svp    = zeros(1,nen);                                         % material point basis functions
  mpData(mp).dSvp   = zeros(nD,nen);                                        % derivative of the basis functions
  mpData(mp).Fn     = eye(3);                                               % previous deformation gradient
  mpData(mp).F      = eye(3);                                               % deformation gradient
  mpData(mp).sig    = zeros(6,1);                                           % Cauchy stress
  mpData(mp).epsEn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
  mpData(mp).epsE   = zeros(6,1);                                           % elastic strain (logarithmic)
  mpData(mp).mCst   = mCst;                                                 % material constants (or internal variables) for constitutive model
  if mpC(mp,1) == max(mpC(:,1)) && mpC(mp,2) == 20 + abs(max(20-mpC(:,2)))
      mpData(mp).fp = [0 P].';
      count = count + 1;
  %   if abs(mpC(mp,1)-(l-lp(mp,1)))<1e-3/a && abs(mpC(mp,2)-(ly-d))<lp(mp,2) %(ly-d))<(lp(mp,2)+1e-3/a)
%     mpData(mp).fp   = [0 P/2].';                                            % point forces at material points (end load)  
  else
    mpData(mp).fp   = zeros(nD,1);                                          % point forces at material points
  end
  mpData(mp).u      = zeros(nD,1);                                          % material point displacements
  if mpData(mp).mpType == 2
    mpData(mp).lp     = lp(mp,:);                                           % material point domain lengths (GIMP)
    mpData(mp).lp0    = lp(mp,:);                                           % initial material point domain lengths (GIMP)
    mpData(mp).R = zeros(3);
  else
    mpData(mp).lp     = zeros(1,nD);                                        % material point domain lengths (MPM)
    mpData(mp).lp0    = zeros(1,nD);                                        % initial material point domain lengths (MPM)
  end
  mpData(mp).mpSte = mpSte;
end
end