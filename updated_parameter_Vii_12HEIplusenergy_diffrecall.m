function viinew = updated_parameter_Vii_12HEIplusenergy_diffrecall(sum_weightrecall,i,r,theta,V,qqa,qqb,qqaw,qqbw,n);
%
% In order to enhance computation efficiency, compute
% qqa, qqb, qqaw, qqbw in the main function instead of computing them inside the
% function.
% 
% call function: formGofSigmaeV_12HEIplusenergy
%
% INPUT:
%      sum_weightrecall: sum of the product of survey weights and # of recalls
%      i:                diagonal position
%      r:                current value of r  
%      theta:            current value of theta
%      V:                matrix V 
%      qqa:              (W-XB-U) for the 1st 24h recall
%      qqb:              (W-XB-U) for the 2nd 24h recall    
%      qqaw:             qqa times weight
%      qqbw:             qqb times weight
%      n:                number of individuals
%
% OUTPUT:
%      Do the Metropolis Step for (i,i)th entry of matrix V
%
Vcurr         = V;
Vcand         = V;
viicurr       = V(i,i);
viicand       = V(i,i) + (0.4 .* (rand(1,1) - 0.5));
Vcand(i,i)    = viicand;
GofSigmaecurr = formGofSigmaeV_12HEIplusenergy(r,theta,Vcurr,qqa,qqb,qqaw,qqbw);
GofSigmaecand = formGofSigmaeV_12HEIplusenergy(r,theta,Vcand,qqa,qqb,qqaw,qqbw);
gg            = GofSigmaecand - GofSigmaecurr;
gg            = gg - (0.5*sum_weightrecall.* log(viicand^2)) + (0.5*sum_weightrecall.* log(viicurr^2));
gg = min([1 exp(gg)*((viicand>=-3)&(viicand<=3))]);
vv            = rand(1,1);
viinew       = (viicand .* (vv < gg)) + (viicurr .* (vv > gg));
