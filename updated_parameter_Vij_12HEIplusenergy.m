function vijnew = updated_parameter_Vij_12HEIplusenergy(i,j,r,theta,V,qqa,qqb,qqaw,qqbw,n);
%
% In order to enhance computation efficiency, compute
% qqa, qqb, qqaw, qqbw in the main function instead of computing them inside the
% function.
% 
% call function: formGofSigmaeV_12HEIplusenergy
%
% INPUT:
%      i:     row position
%      j:     column position
%      r:     current value of r  
%      theta: current value of theta
%      V:     matrix V 
%      qqa:   (W-XB-U) for the 1st 24h recall
%      qqb:   (W-XB-U) for the 2nd 24h recall    
%      qqaw:  qqa times weight
%      qqbw:  qqb times weight
%      n:     number of individuals
%
% OUTPUT:
%      Do the Metropolis Step for (i,j)th entry of matrix V
%
Vcurr         = V;
Vcand         = V;
vijcurr       = V(i,j);
vijcand       = V(i,j) + (0.4 .* (rand(1,1) - 0.5));
Vcand(i,j)    = vijcand;
GofSigmaecurr = formGofSigmaeV_12HEIplusenergy(r,theta,Vcurr,qqa,qqb,qqaw,qqbw);
GofSigmaecand = formGofSigmaeV_12HEIplusenergy(r,theta,Vcand,qqa,qqb,qqaw,qqbw);
gg            = GofSigmaecand - GofSigmaecurr;
gg = min([1 exp(gg)*((vijcand>=-3)&(vijcand<=3))]);
vv            = rand(1,1);
vijnew       = (vijcand .* (vv < gg)) + (vijcurr .* (vv > gg));

