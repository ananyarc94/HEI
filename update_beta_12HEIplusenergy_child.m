function betanew = update_beta_12HEIplusenergy_child(recall,weight,iSigmau,iSigmae,prior_beta_cov, ...
                             Xtildei1,Xtildei2,Utildei,Wtildei,n,beta,prior_beta_mean);
%
% The method of construction is the following. 
% W_{ijk} denotes the ith person, the jth VARIABLE, and the kth REPLICATE. 
% On the other hand, Xtilde(ii,:,jj) is the design matrix for the jth
% variable. Similarly, Utildei(:,jj) is the random effects for the jth
% variable. 
%
% INPUT:
%      recall:          number of 24h recalls for each individual
%      weight:          survey weights  
%      iSigmau:         inverse of Sigmau
%      iSigmae:         inverse of Sigmae
%      prior_beta_cov:  prior covariance matrix of beta
%      Xtildei1:        design matrix for the 1st 24h recall
%      Xtildei2:        design matrix for the 2nd 24h recall
%      Utildei:         current Utildei
%      Wtildei:         current Wtildei
%      n:               number of individuals
%      beta:            current beta  
%      prior_beta_mean: prior mean of beta
% 
% OUTPUT:
%      updated beta
%
betanew = 0 .* beta;
mm      = size(beta,1);

for jj = 1:19; % Variable number
    C1 =  (inv(prior_beta_cov(:,:,jj)) * prior_beta_mean(:,jj)) + ...
                          (iSigmae(jj,jj) .* (((weight*ones(1,mm)) .* Xtildei1(:,:,jj))'  * (Wtildei(:,jj,1) - Utildei(:,jj)))) + ...
                          (iSigmae(jj,jj) .* ((((recall==2) .* weight*ones(1,mm)) .* Xtildei2(:,:,jj))'  * (Wtildei(:,jj,2) - Utildei(:,jj))));
           for ll = 1:19;
            if abs(ll - jj) > 0;
                qq1 = (Wtildei(:,ll,1) - (Xtildei1(:,:,ll) * beta(:,ll)) - Utildei(:,ll)) ;                      
                qq2 = (Wtildei(:,ll,2) - (Xtildei2(:,:,ll) * beta(:,ll)) - Utildei(:,ll)) .*(recall==2);                
                C1 = C1 + (iSigmae(jj,ll) .*( (((weight*ones(1,mm)) .* Xtildei1(:,:,jj))' * qq1)+ ...
                                              (((weight*ones(1,mm)) .* Xtildei2(:,:,jj))' * qq2) ) );
            end;
        end;    
        
   C2 = inv(inv(prior_beta_cov(:,:,jj)) + (iSigmae(jj,jj) .* ( (( weight*ones(1,mm)) .* Xtildei1(:,:,jj))'  * Xtildei1(:,:,jj)+...
                                                               (((recall==2).* weight*ones(1,mm)) .* Xtildei2(:,:,jj))'  * Xtildei2(:,:,jj) ) ));       
   C2 = (C2 + C2')./2;
    betanew(:,jj) = (C2 * C1) + (sqrtm(C2) * randn(mm,1));
end;