function Wtildeinew = gen_Wtildei_12HEIplusenergy_child(recall,Wtildei,beta,Xtildei1,Xtildei2,Utildei,n,iSigmae,Wistar,numgen);
%
% The method of construction is the following. 
% W_{ijk} denotes the ith person, the jth VARIABLE, and the kth REPLICATE. 
% On the other hand, Xtilde(ii,:,jj) is the design matrix for the jth
% variable. Similarly, Utildei(:,jj) is the random effects for the jth
% variable. 
%
% call function: gen_truncated_normals
%
% INPUT:
%      recall:   number of 24h recalls for each individual
%      Wtildei:  current Wtildei
%      beta:     current beta  
%      Xtildei1: design matrix for the 1st 24h recall
%      Xtildei2: design matrix for the 2nd 24h recall
%      Utildei:  current Utildei
%      n:        number of individuals
%      iSigmae:  inverse of Sigmae
%      Wistar:   indicators of consumption of dietary components
%      numgen:   number of times you try to get truncated normals(recommended = 50)
%
% OUTPUT:
%      updated Wtildei
%
Wtildeinew  = Wtildei;   
Wtildeio    = Wtildei;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update W_{i1k} W_{i3k} W_{i5k} W_{i7k} W_{i9k} W_{i11k}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for varnum = 1:2:11; % The variable you want to generate from the truncated normal
kk  = 1;   
    C2      = 1 ./ iSigmae(varnum,varnum);
    C1      = iSigmae(varnum,varnum) .* ((Xtildei1(:,:,varnum) * beta(:,varnum)) + Utildei(:,varnum));
    for jj  = 1:size(Xtildei1,3);
        if abs(jj - varnum) > 0;
            qq = (Wtildei(:,jj,kk) - (Xtildei1(:,:,jj) * beta(:,jj)) - Utildei(:,jj));
            C1 = C1 - (iSigmae(varnum,jj) .* qq);
        end;
    end;
    mu      = C2 .* C1;
    sigma   = sqrt(C2);
    startxi = mu./sigma;
    genww1  = gen_truncated_normals(-mu./sigma,-startxi,numgen);
    genww2  = gen_truncated_normals(mu./sigma,-startxi,numgen);
    if kk  == 1;
        Wtildeinew(:,varnum,kk) = mu + (sigma .* ((Wistar(:,varnum) .* genww1) - ((1 - Wistar(:,varnum)) .* genww2)));
    end;
    if kk  == 2;
        Wtildeinew(:,varnum,kk) = (mu + (sigma .* ((Wistar(:,varnum+1) .* genww1) - ((1 - Wistar(:,varnum+1)) .* genww2)))).* (recall ==2);
    end;
kk  = 2;   
    C2      = 1 ./ iSigmae(varnum,varnum);
    C1      = iSigmae(varnum,varnum) .* ((Xtildei2(:,:,varnum) * beta(:,varnum)) + Utildei(:,varnum));
    for jj  = 1:size(Xtildei2,3);
        if abs(jj - varnum) > 0;
            qq = (Wtildei(:,jj,kk) - (Xtildei2(:,:,jj) * beta(:,jj)) - Utildei(:,jj));
            C1 = C1 - (iSigmae(varnum,jj) .* qq);
        end;
    end;
    mu      = C2 .* C1;
    sigma   = sqrt(C2);
    startxi = mu./sigma;
    genww1  = gen_truncated_normals(-mu./sigma,-startxi,numgen);
    genww2  = gen_truncated_normals(mu./sigma,-startxi,numgen);
    if kk  == 1;
        Wtildeinew(:,varnum,kk) = mu + (sigma .* ((Wistar(:,varnum) .* genww1) - ((1 - Wistar(:,varnum)) .* genww2)));
    end;
    if kk  == 2;
        Wtildeinew(:,varnum,kk) = (mu + (sigma .* ((Wistar(:,varnum+1) .* genww1) - ((1 - Wistar(:,varnum+1)) .* genww2)))).* (recall ==2);
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update W_{ilk},l=2,4,6,8,10,12. You only do this for the cases that W_{ilk} < 0, l= 1,3,5,7,9,11;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for varnum = 2:2:12;  % The variable you want to generate from the truncated normal
kk=1;
    C2      = 1 ./ iSigmae(varnum,varnum);
    C1      = iSigmae(varnum,varnum) .* ((Xtildei1(:,:,varnum) * beta(:,varnum)) + Utildei(:,varnum));
    for jj  = 1:size(Xtildei1,3);
        if abs(jj - varnum) > 0;
            qq = (Wtildeinew(:,jj,kk) - (Xtildei1(:,:,jj) * beta(:,jj)) - Utildei(:,jj));
            C1 = C1 - (iSigmae(varnum,jj) .* qq);
        end;
    end;
    mu      = C2 .* C1;
    sigma   = sqrt(C2);
    Wtildeinew(:,varnum,kk) = mu + (sigma .* randn(n,1));
kk=2;
    C2      = 1 ./ iSigmae(varnum,varnum);
    C1      = iSigmae(varnum,varnum) .* ((Xtildei2(:,:,varnum) * beta(:,varnum)) + Utildei(:,varnum));
    for jj  = 1:size(Xtildei2,3);
        if abs(jj - varnum) > 0;
            qq = (Wtildeinew(:,jj,kk) - (Xtildei2(:,:,jj) * beta(:,jj)) - Utildei(:,jj));
            C1 = C1 - (iSigmae(varnum,jj) .* qq);
        end;
    end;
    mu      = C2 .* C1;
    sigma   = sqrt(C2);
    Wtildeinew(:,varnum,kk) = mu + (sigma .* randn(n,1));

 % Now make sure the actual observed data do not change
 for kk = 1:2;
    if kk == 1;
        Wtildeinew(:,varnum,kk) = (   Wtildei(:,varnum,kk) .* Wistar(:,varnum-1)) + ...
                               (Wtildeinew(:,varnum,kk) .* (1-Wistar(:,varnum-1)));
    end;
    if kk == 2;
        Wtildeinew(:,varnum,kk) = ((   Wtildei(:,varnum,kk) .* Wistar(:,varnum)) + ...
                               (Wtildeinew(:,varnum,kk) .* (1-Wistar(:,varnum)))).* (recall ==2);
    end;
 end;
end;