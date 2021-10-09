function thetainew = updated_parameter_theta_12HEIplusenergy(i,r,theta,V,qqa,qqb,qqaw,qqbw,n);
%
% In order to enhance computation efficiency, compute
% qqa, qqb, qqaw, qqbw in the main function instead of computing them inside the
% function.
% 
% call function: formGofSigmaeV_12HEIplusenergy
%
% INPUT:
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
%      Do the Metropolis Step for theta_i
%
thetapossible = pi .* linspace(-0.99,0.99,41)'; %'
thetamin      = min([thetapossible' theta(i,1)]);
thetamax      = max([thetapossible' theta(i,1)]);
spacing       = thetapossible(2,1) - thetapossible(1,1);
thetacurr  = theta;
thetacand  = theta;
thetaicurr     = theta(i,1);

if thetaicurr <= thetamin;
    ss        = rand(1,1);
    thetaicand = (thetaicurr .* (ss <= 0.33)) + ((thetaicurr + spacing) .* (ss > 0.33) .* (ss <= 0.66)) + ...
                                            + ((thetaicurr + (2 .* spacing)) .* (ss > 0.66));
end;
if thetaicurr >= thetamax;
    ss        = rand(1,1);
    thetaicand = (thetaicurr .* (ss <= 0.33)) + ((thetaicurr - spacing) .* (ss > 0.33) .* (ss <= 0.66)) + ...
                                    + ((thetaicurr - (2 .* spacing)) .* (ss > 0.66));
end;
if thetaicurr  > thetamin;
    if thetaicurr < thetamax;
    ss        = rand(1,1);
    thetaicand = (thetaicurr .* (ss <= 0.33)) + ((thetaicurr + spacing) .* (ss > 0.33) .* (ss <= 0.66)) + ...
                                    + ((thetaicurr - spacing) .* (ss > 0.66));
    end;
end;
thetacand(i,1) = thetaicand;

GofSigmaecurr = formGofSigmaeV_12HEIplusenergy(r,thetacurr,V,qqa,qqb,qqaw,qqbw);
GofSigmaecand = formGofSigmaeV_12HEIplusenergy(r,thetacand,V,qqa,qqb,qqaw,qqbw);
gg            = min([1 exp(GofSigmaecand - GofSigmaecurr)]);
ss            = rand(1,1);
thetainew      = (thetaicand .* (ss <= gg)) + (thetaicurr .* (ss > gg));