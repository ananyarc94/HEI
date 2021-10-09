function rinew = updated_parameter_r_12HEIplusenergy_diffrecall(sum_weightrecall,i,r,theta,V,qqa,qqb,qqaw,qqbw,n);
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
%      Do the Metropolis Step for r_i
%
rpossible = linspace(-0.99,0.99,41)'; %'
rmin      = min([rpossible' r(i,1)]);
rmax      = max([rpossible' r(i,1)]);
spacing   = 0.0495;  %rpossible(2,1) - rpossible(1,1);
rcurr = r;
rcand = r;
ricurr = r(i,1);

if ricurr <= rmin;
    ss    = rand(1,1);
    ricand = (ricurr .* (ss <= 0.33)) + ((ricurr + spacing) .* (ss > 0.33) .* (ss <= 0.66)) + ...
                                    + ((ricurr + (2 .* spacing)) .* (ss > 0.66));
end;
if ricurr >= rmax;
    ss    = rand(1,1);
    ricand = (ricurr .* (ss <= 0.33)) + ((ricurr - spacing) .* (ss > 0.33) .* (ss <= 0.66)) + ...
                                    + ((ricurr - (2 .* spacing)) .* (ss > 0.66));
end;
if ricurr > rmin;
    if ricurr < rmax;
    ss    = rand(1,1);
    ricand = (ricurr .* (ss <= 0.33)) + ((ricurr + spacing) .* (ss > 0.33) .* (ss <= 0.66)) + ...
                                    + ((ricurr - spacing) .* (ss > 0.66));
    end;
end;
rcand(i,1) = ricand;

GofSigmaecurr = formGofSigmaeV_12HEIplusenergy(rcurr,theta,V,qqa,qqb,qqaw,qqbw);
GofSigmaecand = formGofSigmaeV_12HEIplusenergy(rcand,theta,V,qqa,qqb,qqaw,qqbw);
gg            = GofSigmaecand - GofSigmaecurr;
gg            = gg - (0.5*sum_weightrecall .* log(1 - (ricand .^ 2))) + (0.5*sum_weightrecall  .* log(1 - (ricurr .^ 2)));
gg = min([1 exp(gg)*((ricand>=-1)&(ricand<=1))]);
ss            = rand(1,1);
rinew          = (ricand .* (ss <= gg)) + (ricurr .* (ss > gg));