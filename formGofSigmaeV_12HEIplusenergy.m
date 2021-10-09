function GofSigmaeV = formGofSigmaeV_12HEIplusenergy(r,theta,V,qqa,qqb,qqaw,qqbw);
%
% This calculates the part inside the exponential of function g() in the paper 
%
% In order to enhance computation efficiency, compute
% qqa, qqb, qqaw, qqbw in the main function instead of computing them inside the
% function.
% 
% INPUT:
%      r:                current value of r  
%      theta:            current value of theta
%      V:                matrix V 
%      qqa:              (W-XB-U) for the 1st 24h recall
%      qqb:              (W-XB-U) for the 2nd 24h recall    
%      qqaw:             qqa times weight
%      qqbw:             qqb times weight
%
% OUTPUT:
%      the part inside the exponential of function g() in the paper
%
V(1,1) = 1;
V(2,1) = 0;

V(3,1) = r(1,1) .* sin(theta(1,1));
V(3,2) = r(1,1) .* cos(theta(1,1));
 
for k = 2:5;

 V(2*k+1,1) = r(k,1) .* sin(theta((k-1)^2+1,1));
 
  for j = 2:(2*k-1);
    V(2*k+1,j) = r(k,1) .* sin(theta((k-1)^2+j,1)) ;
     for jj = 1: (j-1);
      V(2*k+1,j) = V(2*k+1,j).* cos(theta((k-1)^2+jj,1));
     end;
  end;
 
 V(2*k+1,2*k) = r(k,1);
 for jj = 1: (2*k-1);
  V(2*k+1,2*k) = V(2*k+1,2*k).* cos(theta((k-1)^2+jj,1));
 end;

end;

V(3,3)  = sqrt(1-r(1,1)^2);
V(5,5)  = sqrt(1-r(2,1)^2);
V(7,7)  = sqrt(1-r(3,1)^2);
V(9,9)  = sqrt(1-r(4,1)^2);
V(11,11)= sqrt(1-r(5,1)^2);

for k = 3:2:11
V(k+1,k) = -(sum(V(k,1:(k-1)).*V(k+1,1:(k-1))))./ V(k,k);
end;

Sigmae = V*V';
iSigmae      = inv(Sigmae);

tempMat = qqaw' * qqa + qqbw'*qqb;
tempMat = iSigmae.*tempMat;
GofSigmaeV = -0.5*(sum(sum(tempMat)));