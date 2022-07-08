function [stiffness]=assemble(stiffness,k,index)
 edof = length(index);
 for i=1:edof
     ii=index(i);
     for j=1:edof
         jj=index(j);
         stiffness(ii,jj)=stiffness(ii,jj)+k(i,j);
     end
 end