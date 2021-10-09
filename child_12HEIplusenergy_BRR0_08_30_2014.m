%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ( unsetenv DISPLAY ; matlab < child_12HEIplusenergy_BRR0.m >& /dev/null ) &
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the MCMC analysis on fake data for 12 HEI in the two-part model
%
% Call in:
% formGofSigmaeV_12HEIplusenergy.m
% gen_truncated_normals.m
% gen_Wtildei_12HEIplusenergy_child.m
% update_beta_12HEIplusenergy_child.m
% updated_parameter_r_12HEIplusenergy_diffrecall.m
% updated_parameter_theta_12HEIplusenergy.m
% updated_parameter_Vii_12HEIplusenergy_diffrecall.m
% updated_parameter_Vij_12HEIplusenergy.m
%
% INPUT:
%      NHANES_noname_-1.csv
%
% OUTPUT:
%     child_Sigmau_mean.mat 
%     child_Sigmae_mean.mat 
%     child_beta_mean.mat      
%     child_Sigmau.mat 
%     child_Sigmae.mat
%     child_beta.mat   
%     q1.mat
%     q2.mat                 % Table 2 in the paper
%     q3.mat                 % Table 3 in the paper
%     adjustedintake_R.txt   % Table 4 in the paper
%     score_R_model.mat      % part of Table 5 in the paper
%     q4_lt50.mat            % part of Table 6 in the paper
%     q4_gt50.mat            % part of Table 6 in the paper 
%     Figure 1,2 and 3 in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%set to your own path
path(path,['C:\Users\anany\OneDrive\Desktop\Research\HEI\NHANES_Data_Code_AOAS'])
%You may select another starting seed
myseed = 6309021;
format compact;
rand('state',myseed);
randn('state',myseed);

tic
nMCMC     = 700;
nburn     = 200;

thecount2 = 0;

diary('Test_Output/child_12HEIplusenergy_BRR0.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set the parameters');
data = csvread('NHANES_noname_-1.csv');                               
                
child           = ((data(:,8)>=2) & (data(:,8)<=8));          
data            = data(child,:);                                       
                                                                      
yy              = data(~ismember(data(:,1),data(data(:,52)==2,1)),:);  
zz              = data(ismember(data(:,1),data(data(:,52)==2,1)),:);   

n1 = size(yy,1); 
n2 = size(zz,1)/2;
n  = n1+n2;

recall = [(2*ones(n2,1))', ones(n1,1)']';      % number of recalls per person
DAY =zz(:,52);
DAY1 = [zz(DAY==1,52)',yy(:,52)']';  
DAY2 = [zz(DAY==2,52)',yy(:,52)']';  
day1 = (DAY1 ==2);  
day2 = (DAY2 ==2);  

WKEND1 = [zz(DAY==1,54)',yy(:,54)']';  
WKEND2 = [zz(DAY==2,54)',yy(:,54)']';  
wkend1 = (WKEND1 ==1);  
wkend2 = (WKEND2 ==1);  
age	             = [zz(DAY==2,8)',yy(:,8)']';            
sex	             = [zz(DAY==2,6)',yy(:,6)']';            
sex	             = (sex ==1); 
race	             = [zz(DAY==2,11)',yy(:,11)']';                                           
race2              = (race ==2); 

weight             = [zz(DAY==2,17)',yy(:,17)']';

HEI1_F_TOT	       = [zz(DAY==1,173)',yy(:,173)']';
HEI2_F_TOT	       = [zz(DAY==2,173)',zeros(n1,1)']';
HEI1_F_WHL	       = [zz(DAY==1,174)',yy(:,174)']';
HEI2_F_WHL	       = [zz(DAY==2,174)',zeros(n1,1)']';
HEI1_F_other            = HEI1_F_TOT - HEI1_F_WHL;              
HEI2_F_other            = HEI2_F_TOT -HEI2_F_WHL;
HEI1_V_TOT	       = [zz(DAY==1,176)',yy(:,176)']';
HEI2_V_TOT	       = [zz(DAY==2,176)',zeros(n1,1)']';
HEI1_V_DOL	       = [zz(DAY==1,177)',yy(:,177)']';
HEI2_V_DOL	       = [zz(DAY==2,177)',zeros(n1,1)']';
HEI1_V_other            = HEI1_V_TOT - HEI1_V_DOL;              
HEI2_V_other            = HEI2_V_TOT - HEI2_V_DOL;
HEI1_G_WHL	       = [zz(DAY==1,179)',yy(:,179)']';
HEI2_G_WHL	       = [zz(DAY==2,179)',zeros(n1,1)']';
HEI1_MILK	       = [zz(DAY==1,180)',yy(:,180)']';
HEI2_MILK	       = [zz(DAY==2,180)',zeros(n1,1)']';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if there are only a few zeroes for a component that is to be treated
%like a non-episodically-consumed food, our standard practice is to replace all
%zeroes with 1/2 the minimum value of the nonzeros for that food.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HEI1_G_TOT	              = [zz(DAY==1,178)',yy(:,178)']';
HEI1_G_TOT(HEI1_G_TOT==0) = 0.5 .* min(HEI1_G_TOT(HEI1_G_TOT>0));
temp          = zz(DAY==2,178);
temp(temp==0) = 0.5 .* min(temp(temp>0));
HEI2_G_TOT	  = [temp',zeros(n1,1)']';

%%%%%%%%%%%%%%%%%%%%%%%%%%
HEI1_G_other            = HEI1_G_TOT - HEI1_G_WHL;  
HEI1_G_other(HEI1_G_other==0) = 0.5 .* min(HEI1_G_other(HEI1_G_other>0));            
HEI2_G_other            = HEI2_G_TOT - HEI2_G_WHL;
HEI2_G_other(HEI2_G_other==0) = 0.5 .* min(HEI2_G_other(HEI2_G_other>0));
HEI2_G_other	  = [HEI2_G_other(1:n2,1)',zeros(n1,1)']';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HEI1_M_TOT	              = [zz(DAY==1,181)',yy(:,181)']';
HEI1_M_TOT(HEI1_M_TOT==0) = 0.5 .* min(HEI1_M_TOT(HEI1_M_TOT>0));
temp          = zz(DAY==2,181);
temp(temp==0) = 0.5 .* min(temp(temp>0));
HEI2_M_TOT	  = [temp', zeros(n1,1)']';

HEI1_OIL	          = [zz(DAY==1,182)',yy(:,182)']';
HEI1_OIL(HEI1_OIL==0) = 0.5 .* min(HEI1_OIL(HEI1_OIL>0));
temp          = zz(DAY==2,182);
temp(temp==0) = 0.5 .* min(temp(temp>0));
HEI2_OIL	  = [temp', zeros(n1,1)']';


HEI1_SFAT	            = [zz(DAY==1,183)',yy(:,183)']';
HEI1_SFAT(HEI1_SFAT==0) = 0.5 .* min(HEI1_SFAT(HEI1_SFAT>0));
temp          = zz(DAY==2,183);
temp(temp==0) = 0.5 .* min(temp(temp>0));
HEI2_SFAT	  = [temp', zeros(n1,1)']';

HEI1_SODI	            = [zz(DAY==1,184)',yy(:,184)']';
HEI1_SODI(HEI1_SODI==0) = 0.5 .* min(HEI1_SODI(HEI1_SODI>0));
temp          = zz(DAY==2,184);
temp(temp==0) = 0.5 .* min(temp(temp>0));
HEI2_SODI	  = [temp', zeros(n1,1)']';


HEI1_SOFAAS_C	              = [zz(DAY==1,185)',yy(:,185)']';
HEI1_SOFAAS_C(HEI1_SOFAAS_C==0) = 0.5 .* min(HEI1_SOFAAS_C(HEI1_SOFAAS_C>0));
temp          = zz(DAY==2,185);
temp(temp==0) = 0.5 .* min(temp(temp>0));
HEI2_SOFAAS_C = [temp', zeros(n1,1)']';

energy1            = [zz(DAY==1,55)', yy(:,55)']';
energy2            = [zz(DAY==2,55)', zeros(n1,1)']';

%%%%%%%%%%%%%%%%%%%%%%%%%%
energy1_other            = energy1 - 9.*HEI1_SFAT-HEI1_SOFAAS_C;  
energy1_other(energy1_other==0) = 0.5 .* min(energy1_other(energy1_other>0));            
energy2_other            = energy2 - 9.*HEI2_SFAT-HEI2_SOFAAS_C;
energy2_other(energy2_other==0) = 0.5 .* min(energy2_other(energy2_other>0));
energy2_other	  = [energy2_other(1:n2,1)',zeros(n1,1)']'; %'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform, delete outliers and standardize the recalls for 6 epi-foods , 
% make sure zeros stay as zeros
% Transform and standardize the recalls for 3 non-epi ,3 nutr, energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 lam = zeros(13,1);
 lam(1,1) =  0.38;
 lam(2,1) =  0.42;
 lam(3,1) =  0.38;
 lam(4,1) =  0.18;
 lam(5,1) =  0.32;
 lam(6,1) =  0.53;
 lam(7,1) =  0.35;
 lam(8,1) =  0.56;
 lam(9,1) =  0.38;
 lam(10,1)=  0.32;
 lam(11,1)=  0.28;
 lam(12,1)=  0.28;
 lam(13,1)=  0.36;

Wistar           = zeros(n,24);
Wistar(:,1)      = (HEI1_F_other  > 0);	       
Wistar(:,2)      = (HEI2_F_other  > 0);	       
Wistar(:,3)      = (HEI1_F_WHL    > 0);
Wistar(:,4)      = (HEI2_F_WHL    > 0);	       
Wistar(:,5)      = (HEI1_V_other  > 0);	       
Wistar(:,6)      = (HEI2_V_other  > 0);	       
Wistar(:,7)      = (HEI1_V_DOL    > 0);       
Wistar(:,8)      = (HEI2_V_DOL    > 0);
Wistar(:,9)      = (HEI1_G_WHL    > 0);
Wistar(:,10)     = (HEI2_G_WHL    > 0);       
Wistar(:,11)     = (HEI1_MILK	    > 0);       
Wistar(:,12)     = (HEI2_MILK     > 0);

Wistar(:,13)     = (HEI1_G_other  > 0);	       
Wistar(:,14)     = (HEI2_G_other  > 0);	       
Wistar(:,15)     = (HEI1_M_TOT    > 0);
Wistar(:,16)     = (HEI2_M_TOT    > 0);	       
Wistar(:,17)     = (HEI1_OIL      > 0);	       
Wistar(:,18)     = (HEI2_OIL      > 0);       
Wistar(:,19)     = (HEI1_SFAT     > 0);       
Wistar(:,20)     = (HEI2_SFAT     > 0);
Wistar(:,21)     = (HEI1_SODI     > 0);
Wistar(:,22)     = (HEI2_SODI     > 0);       
Wistar(:,23)     = (HEI1_SOFAAS_C > 0);       
Wistar(:,24)     = (HEI2_SOFAAS_C > 0);

Wi021        = (HEI1_F_other .^lam(1,1)-1)/lam(1,1);
Wi022        = (HEI2_F_other .^lam(1,1)-1)/lam(1,1);
Wi041        = (HEI1_F_WHL   .^lam(2,1)-1)/lam(2,1);
Wi042        = (HEI2_F_WHL   .^lam(2,1)-1)/lam(2,1);
Wi061        = (HEI1_V_other .^lam(3,1)-1)/lam(3,1);
Wi062        = (HEI2_V_other .^lam(3,1)-1)/lam(3,1);
Wi081        = (HEI1_V_DOL   .^lam(4,1)-1)/lam(4,1);
Wi082        = (HEI2_V_DOL   .^lam(4,1)-1)/lam(4,1);
Wi101        = (HEI1_G_WHL   .^lam(5,1)-1)/lam(5,1);
Wi102        = (HEI2_G_WHL   .^lam(5,1)-1)/lam(5,1);
Wi121        = (HEI1_MILK    .^lam(6,1)-1)/lam(6,1);
Wi122        = (HEI2_MILK    .^lam(6,1)-1)/lam(6,1);

Wi131        = (HEI1_G_other .^lam(7,1)-1)/lam(7,1);
Wi132        = (HEI2_G_other .^lam(7,1)-1)/lam(7,1);
Wi141        = (HEI1_M_TOT   .^lam(8,1)-1)/lam(8,1);
Wi142        = (HEI2_M_TOT   .^lam(8,1)-1)/lam(8,1);
Wi151        = (HEI1_OIL     .^lam(9,1)-1)/lam(9,1);
Wi152        = (HEI2_OIL     .^lam(9,1)-1)/lam(9,1);
Wi161        = (HEI1_SFAT    .^lam(10,1)-1)/lam(10,1);
Wi162        = (HEI2_SFAT    .^lam(10,1)-1)/lam(10,1);
Wi171        = (HEI1_SODI    .^lam(11,1)-1)/lam(11,1);
Wi172        = (HEI2_SODI    .^lam(11,1)-1)/lam(11,1);
Wi181        = (HEI1_SOFAAS_C.^lam(12,1)-1)/lam(12,1);
Wi182        = (HEI2_SOFAAS_C.^lam(12,1)-1)/lam(12,1);

Wi191        = (energy1_other .^lam(13,1)-1)./lam(13,1);
Wi192        = (energy2_other .^lam(13,1)-1)./lam(13,1);

q1=prctile([Wi021(HEI1_F_other> 0)' Wi022(HEI2_F_other> 0)']', 25);
q3=prctile([Wi021(HEI1_F_other> 0)' Wi022(HEI2_F_other> 0)']', 75);
iqr = q3-q1;
index021 = (HEI1_F_other> 0)&(Wi021>q3+3*iqr );
index022 = (HEI2_F_other> 0)&(Wi022>q3+3*iqr );

q1=prctile([Wi041(HEI1_F_WHL> 0)' Wi042(HEI2_F_WHL> 0)']', 25);
q3=prctile([Wi041(HEI1_F_WHL> 0)' Wi042(HEI2_F_WHL> 0)']', 75);
iqr = q3-q1;
index041 = (HEI1_F_WHL> 0)&(Wi041>q3+3*iqr) ;
index042 = (HEI2_F_WHL> 0)&(Wi042>q3+3*iqr) ;

q1=prctile([Wi061(HEI1_V_other> 0)' Wi062(HEI2_V_other> 0)']', 25);
q3=prctile([Wi061(HEI1_V_other> 0)' Wi062(HEI2_V_other> 0)']', 75);
iqr = q3-q1;
index061 = (HEI1_V_other> 0)&(Wi061>q3+3*iqr) ;
index062 = (HEI2_V_other> 0)&(Wi062>q3+3*iqr) ;

q1=prctile([Wi081(HEI1_V_DOL> 0)' Wi082(HEI2_V_DOL> 0)']', 25);
q3=prctile([Wi081(HEI1_V_DOL> 0)' Wi082(HEI2_V_DOL> 0)']', 75);
iqr = q3-q1;
index081 = (HEI1_V_DOL> 0)&(Wi081>q3+3*iqr) ;
index082 = (HEI2_V_DOL> 0)&(Wi082>q3+3*iqr) ;

q1=prctile([Wi101(HEI1_G_WHL> 0)' Wi102(HEI2_G_WHL> 0)']', 25);
q3=prctile([Wi101(HEI1_G_WHL> 0)' Wi102(HEI2_G_WHL> 0)']', 75);
iqr = q3-q1;
index101 = (HEI1_G_WHL> 0)&(Wi101>q3+3*iqr) ;
index102 = (HEI2_G_WHL> 0)&(Wi102>q3+3*iqr) ;

q1=prctile([Wi121(HEI1_MILK> 0)' Wi122(HEI2_MILK> 0)']', 25);
q3=prctile([Wi121(HEI1_MILK> 0)' Wi122(HEI2_MILK> 0)']', 75);
iqr = q3-q1;
index121 = (HEI1_MILK> 0)&(Wi121>q3+3*iqr) ;
index122 = (HEI2_MILK> 0)&(Wi122>q3+3*iqr) ;

q1=prctile([Wi131(HEI1_G_other> 0)' Wi132(HEI2_G_other> 0)']', 25);
q3=prctile([Wi131(HEI1_G_other> 0)' Wi132(HEI2_G_other> 0)']', 75);
iqr = q3-q1;
index131 = (HEI1_G_other> 0)&(Wi131>q3+3*iqr) ;
index132 = (HEI2_G_other> 0)&(Wi132>q3+3*iqr) ;

q1=prctile([Wi141(HEI1_M_TOT> 0)' Wi142(HEI2_M_TOT> 0)']', 25);
q3=prctile([Wi141(HEI1_M_TOT> 0)' Wi142(HEI2_M_TOT> 0)']', 75);
iqr = q3-q1;
index141 = (HEI1_M_TOT> 0)&(Wi141>q3+3*iqr);
index142 = (HEI2_M_TOT> 0)&(Wi142>q3+3*iqr);

q1=prctile([Wi151(HEI1_OIL> 0)' Wi152(HEI2_OIL> 0)']', 25);
q3=prctile([Wi151(HEI1_OIL> 0)' Wi152(HEI2_OIL> 0)']', 75);
iqr = q3-q1;
index151 = (HEI1_OIL> 0)&(Wi151>q3+3*iqr) ;
index152 = (HEI2_OIL> 0)&(Wi152>q3+3*iqr);

q1=prctile([Wi161(HEI1_SFAT> 0)' Wi162(HEI2_SFAT> 0)']', 25);
q3=prctile([Wi161(HEI1_SFAT> 0)' Wi162(HEI2_SFAT> 0)']', 75);
iqr = q3-q1;
index161 = (HEI1_SFAT> 0)&(Wi161>q3+3*iqr) ;
index162 = (HEI2_SFAT> 0)&(Wi162>q3+3*iqr) ;

q1=prctile([Wi171(HEI1_SODI> 0)' Wi172(HEI2_SODI> 0)']', 25);
q3=prctile([Wi171(HEI1_SODI> 0)' Wi172(HEI2_SODI> 0)']', 75);
iqr = q3-q1;
index171 = (HEI1_SODI> 0)&(Wi171>q3+3*iqr);
index172 = (HEI2_SODI> 0)&(Wi172>q3+3*iqr) ;

q1=prctile([Wi181(HEI1_SOFAAS_C> 0)' Wi182(HEI2_SOFAAS_C> 0)']', 25);
q3=prctile([Wi181(HEI1_SOFAAS_C> 0)' Wi182(HEI2_SOFAAS_C> 0)']', 75);
iqr = q3-q1;
index181 = (HEI1_SOFAAS_C> 0)&(Wi181>q3+3*iqr) ;
index182 = (HEI2_SOFAAS_C> 0)&(Wi182>q3+3*iqr) ;

q1=prctile([Wi191(energy1_other> 0)' Wi192(energy2_other> 0)']', 25);
q3=prctile([Wi191(energy1_other> 0)' Wi192(energy2_other> 0)']', 75);
iqr = q3-q1;
index191 = (energy1_other> 0)&(Wi191>q3+3*iqr) ;
index192 = (energy2_other> 0)&(Wi192>q3+3*iqr) ;

outlierindex = ((index021+index022+index041+index042+index061+index062+index081+index082+index101+index102+index121+index122...
                +index131+index132+index141+index142+index151+index152+index161+index162+index171+index172+index181+index182+index191+index192)>0);
nonoutlierindex = ((index021+index022+index041+index042+index061+index062+index081+index082+index101+index102+index121+index122...
                +index131+index132+index141+index142+index151+index152+index161+index162+index171+index172+index181+index182+index191+index192)==0);

mu  = zeros(13,1);
sig = zeros(13,1);

mu(1,1)  = mean([Wi021(HEI1_F_other> 0 & nonoutlierindex)' Wi022(HEI2_F_other> 0 & nonoutlierindex)']');
sig(1,1) =  std([Wi021(HEI1_F_other> 0 & nonoutlierindex)' Wi022(HEI2_F_other> 0 & nonoutlierindex)']');
temp     = (Wistar(:,1).* sqrt(2) .* (Wi021 - mu(1,1)) ./ sig(1,1)); 
Wi021    = temp(nonoutlierindex);
temp     = (Wistar(:,2).* sqrt(2) .* (Wi022 - mu(1,1)) ./ sig(1,1));
Wi022    = temp(nonoutlierindex);

mu(2,1)  = mean([Wi041(HEI1_F_WHL> 0 & nonoutlierindex)' Wi042(HEI2_F_WHL> 0 & nonoutlierindex)']');
sig(2,1) =  std([Wi041(HEI1_F_WHL> 0 & nonoutlierindex)' Wi042(HEI2_F_WHL> 0 & nonoutlierindex)']');
temp     = Wistar(:,3) .* sqrt(2) .* (Wi041 - mu(2,1)) ./ sig(2,1);
Wi041    = temp(nonoutlierindex);
temp     = Wistar(:,4) .* sqrt(2) .* (Wi042 - mu(2,1)) ./ sig(2,1);
Wi042    = temp(nonoutlierindex);

mu(3,1)  = mean([Wi061(HEI1_V_other> 0 & nonoutlierindex)' Wi062(HEI2_V_other> 0 & nonoutlierindex)']');
sig(3,1) =  std([Wi061(HEI1_V_other> 0 & nonoutlierindex)' Wi062(HEI2_V_other> 0 & nonoutlierindex)']');
temp     = (Wistar(:,5).* sqrt(2) .* (Wi061 - mu(3,1)) ./ sig(3,1)); 
Wi061    = temp(nonoutlierindex);
temp     = (Wistar(:,6).* sqrt(2) .* (Wi062 - mu(3,1)) ./ sig(3,1));
Wi062    = temp(nonoutlierindex);

mu(4,1)  = mean([Wi081(HEI1_V_DOL> 0 & nonoutlierindex)'  Wi082(HEI2_V_DOL> 0 & nonoutlierindex)']');
sig(4,1) =  std([Wi081(HEI1_V_DOL> 0 & nonoutlierindex)'  Wi082(HEI2_V_DOL> 0 & nonoutlierindex)']');
temp     = Wistar(:,7) .* sqrt(2) .* (Wi081 - mu(4,1)) ./ sig(4,1);
Wi081    = temp(nonoutlierindex);
temp     = Wistar(:,8) .* sqrt(2) .* (Wi082 - mu(4,1)) ./ sig(4,1);
Wi082    = temp(nonoutlierindex);

mu(5,1)          = mean([Wi101(HEI1_G_WHL> 0 & nonoutlierindex)'  Wi102(HEI2_G_WHL> 0 & nonoutlierindex)']');
sig(5,1)         =  std([Wi101(HEI1_G_WHL> 0 & nonoutlierindex)'  Wi102(HEI2_G_WHL> 0 & nonoutlierindex)']');
temp     =  Wistar(:,9) .* sqrt(2) .* (Wi101 - mu(5,1)  ) ./ sig(5,1)  ;
Wi101    = temp(nonoutlierindex);
temp     = Wistar(:,10) .* sqrt(2) .* (Wi102 - mu(5,1)  ) ./ sig(5,1)  ;
Wi102    = temp(nonoutlierindex);

mu(6,1)         = mean([Wi121(HEI1_MILK> 0 & nonoutlierindex)'  Wi122(HEI2_MILK> 0 & nonoutlierindex)']');
sig(6,1)        =  std([Wi121(HEI1_MILK> 0 & nonoutlierindex)'  Wi122(HEI2_MILK> 0 & nonoutlierindex)']');
temp     = Wistar(:,11).* sqrt(2) .* (Wi121 - mu(6,1)  ) ./ sig(6,1)  ;
Wi121    = temp(nonoutlierindex);
temp     = Wistar(:,12) .* sqrt(2) .* (Wi122 - mu(6,1)  ) ./ sig(6,1)  ;
Wi122    = temp(nonoutlierindex);

mu(7,1)  = mean([Wi131(HEI1_G_other> 0 & nonoutlierindex)' Wi132(HEI2_G_other> 0 & nonoutlierindex)']');
sig(7,1) =  std([Wi131(HEI1_G_other> 0 & nonoutlierindex)' Wi132(HEI2_G_other> 0 & nonoutlierindex)']');
temp     = sqrt(2) .* (Wi131 - mu(7,1)) ./ sig(7,1); 
Wi131    = temp(nonoutlierindex);
temp     = (Wistar(:,14).* sqrt(2) .* (Wi132 - mu(7,1)) ./ sig(7,1));
Wi132    = temp(nonoutlierindex);

mu(8,1)          = mean([Wi141(nonoutlierindex)'  Wi142(HEI2_M_TOT> 0 & nonoutlierindex)']');
sig(8,1)         =  std([Wi141(nonoutlierindex)'  Wi142(HEI2_M_TOT> 0 & nonoutlierindex)']');
temp     =                 sqrt(2) .* (Wi141 - mu(8,1)  ) ./ sig(8,1)  ;
Wi141    = temp(nonoutlierindex);
temp     = Wistar(:,16) .* sqrt(2) .* (Wi142 - mu(8,1)  ) ./ sig(8,1)  ;
Wi142    = temp(nonoutlierindex);

mu(9,1)          = mean([Wi151(nonoutlierindex)'  Wi152(HEI2_OIL> 0 & nonoutlierindex)']');
sig(9,1)         =  std([Wi151(nonoutlierindex)'  Wi152(HEI2_OIL> 0 & nonoutlierindex)']');
temp     =                 sqrt(2) .* (Wi151 - mu(9,1)  ) ./ sig(9,1)  ;
Wi151    = temp(nonoutlierindex);
temp     = Wistar(:,18) .* sqrt(2) .* (Wi152 - mu(9,1)  ) ./ sig(9,1)  ;
Wi152    = temp(nonoutlierindex);

mu(10,1)          = mean([Wi161(nonoutlierindex)'  Wi162(HEI2_SFAT> 0 & nonoutlierindex)']');
sig(10,1)         =  std([Wi161(nonoutlierindex)'  Wi162(HEI2_SFAT> 0 & nonoutlierindex)']');
temp     =                 sqrt(2) .* (Wi161 - mu(10,1)  ) ./ sig(10,1)  ;
Wi161    = temp(nonoutlierindex);
temp     = Wistar(:,20) .* sqrt(2) .* (Wi162 - mu(10,1)  ) ./ sig(10,1)  ;
Wi162    = temp(nonoutlierindex);

mu(11,1)          = mean([Wi171(nonoutlierindex)'  Wi172(HEI2_SODI> 0 & nonoutlierindex)']');
sig(11,1)         =  std([Wi171(nonoutlierindex)'  Wi172(HEI2_SODI> 0 & nonoutlierindex)']');
temp     =                 sqrt(2) .* (Wi171 - mu(11,1)  ) ./ sig(11,1)  ;
Wi171    = temp(nonoutlierindex);
temp     = Wistar(:,22) .* sqrt(2) .* (Wi172 - mu(11,1)  ) ./ sig(11,1)  ;
Wi172    = temp(nonoutlierindex);

mu(12,1)          = mean([Wi181(nonoutlierindex)'  Wi182(HEI2_SOFAAS_C> 0 & nonoutlierindex)']');
sig(12,1)         =  std([Wi181(nonoutlierindex)'  Wi182(HEI2_SOFAAS_C> 0 & nonoutlierindex)']');
temp     =                 sqrt(2) .* (Wi181 - mu(12,1)  ) ./ sig(12,1)  ;
Wi181    = temp(nonoutlierindex);
temp     = Wistar(:,24) .* sqrt(2) .* (Wi182 - mu(12,1)  ) ./ sig(12,1)  ;
Wi182    = temp(nonoutlierindex);

mu(13,1)  = mean([Wi191(nonoutlierindex)' Wi192(energy2_other> 0 & nonoutlierindex)']');
sig(13,1)  =  std([Wi191(nonoutlierindex)' Wi192(energy2_other> 0 & nonoutlierindex)']');
temp     =                sqrt(2) .* (Wi191 - mu(13,1)  ) ./ sig(13,1)  ;
Wi191    = temp(nonoutlierindex);
temp     = (energy2_other> 0).* sqrt(2) .* (Wi192 - mu(13,1)  ) ./ sig(13,1)  ;
Wi192    = temp(nonoutlierindex);
disp(['   From energy, I subtracted ',num2str(mu(13,1) )]);
disp(['           I then divided by ',num2str(sig(13,1)  ./sqrt(2))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the design matrices. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recall= recall (nonoutlierindex);
n1 = sum(recall==1);
n2 = sum(recall==2);
n  = sum(nonoutlierindex);  % or could be calculated in another way: (n1+n2)
day1 = day1(nonoutlierindex);
day2 = day2(nonoutlierindex);
wkend1 = wkend1(nonoutlierindex);
wkend2 = wkend2(nonoutlierindex);
age     = age(nonoutlierindex);
age     = (age - mean(age)) ./ std(age);
sex     = sex(nonoutlierindex);
race2 = race2(nonoutlierindex);
agesex     = age .* sex;
agerace2   = age .* race2;
sexrace2   = sex .* race2;
agesexrace2= age .* sex .* race2;

Xtildei1          = zeros(n,10,19);
Xtildei2          = zeros(n,10,19);

for j = 1:19;
   Xtildei1(:,:,j) = [ones(n,1) day1 wkend1 age sex race2 agesex agerace2 sexrace2 agesexrace2];
end;
for j = 1:19;
   Xtildei2(:,:,j) = [ones(n,1) day2 wkend2 age sex race2 agesex agerace2 sexrace2 agesexrace2];
end;

weight = weight(nonoutlierindex);
weight = weight/ sum(weight) * n;                                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set the parameters');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the prior for beta and starting values for it
%
% Good starting values for the beta parameters and their covariance
% matrices are available by other means. I ran the consumption program and
% the amount program to get these values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prior_beta_mean = zeros(10,19);
prior_beta_cov  = zeros(10,10,19);
for j = 1:19;
 prior_beta_cov(:,:,j) = 100 .* eye(10);
end;
beta_start = prior_beta_mean;
beta       = beta_start;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the prior for Sigmae
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V       = eye(19);
r       = zeros(11,1);
V(1,1)  = 1;
V(2,1)  = 0;
V(2,2)=1;
V(4,4)=1;
theta   = zeros(121,1);
Sigmae = V*V'; %'
iSigmae    = inv(Sigmae);
Sigmae_start = Sigmae;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the prior for Sigmau and starting values for it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prior_Sigmau_mean  = eye(19);
Sigmau             = prior_Sigmau_mean;
for jj = 1:19;
    for kk = 1:19;
        if abs(jj-kk) > 0;
            Sigmau(jj,kk) = 0.50 .* sqrt(Sigmau(jj,jj) .* Sigmau(kk,kk));
        end;
    end;
end;
prior_Sigmau_doff = 1 + 1 + size(Sigmau,1);
prior_Sigmau_mean = Sigmau;
Sigmau    = prior_Sigmau_mean;
iSigmau    = inv(Sigmau);
Sigmau_start = Sigmau;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set starting values for the Utildei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Utildei = randn(n,size(Sigmau,2)) * sqrtm(Sigmau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get starting values for the W_{i1k} and W_{i2k}, the latter 
% when they are not observed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WtildeiS         = zeros(n,19,2);
WtildeiS(:,2,1)  = Wi021;
WtildeiS(:,2,2)  = Wi022;
WtildeiS(:,4,1)  = Wi041;
WtildeiS(:,4,2)  = Wi042;
WtildeiS(:,6,1)  = Wi061;
WtildeiS(:,6,2)  = Wi062;
WtildeiS(:,8,1)  = Wi081;
WtildeiS(:,8,2)  = Wi082;
WtildeiS(:,10,1) = Wi101;
WtildeiS(:,10,2) = Wi102;
WtildeiS(:,12,1) = Wi121;
WtildeiS(:,12,2) = Wi122;
WtildeiS(:,13,1) = Wi131;
WtildeiS(:,13,2) = Wi132;
WtildeiS(:,14,1) = Wi141;
WtildeiS(:,14,2) = Wi142;
WtildeiS(:,15,1) = Wi151;
WtildeiS(:,15,2) = Wi152;
WtildeiS(:,16,1) = Wi161;
WtildeiS(:,16,2) = Wi162;
WtildeiS(:,17,1) = Wi171;
WtildeiS(:,17,2) = Wi172;
WtildeiS(:,18,1) = Wi181;
WtildeiS(:,18,2) = Wi182;
WtildeiS(:,19,1) = Wi191;
WtildeiS(:,19,2) = Wi192;
Wistar           = Wistar(nonoutlierindex,:);
for j =1:2:11;
 WtildeiS(:,j,1)  = abs(Xtildei1(:,:,j) * beta(:,j) + Utildei(:,j) + randn(n,1));
 WtildeiS(:,j,2)  = abs(Xtildei2(:,:,j) * beta(:,j) + Utildei(:,j) + randn(n,1)).* (recall==2);
 WtildeiS(:,j,1)  = (WtildeiS(:,j,1) .* Wistar(:,j)) - (WtildeiS(:,j,1) .* (1 - Wistar(:,j)));
 WtildeiS(:,j,2)  = (WtildeiS(:,j,2) .* Wistar(:,j+1)) - (WtildeiS(:,j,2) .* (1 - Wistar(:,j+1)));
end;

numgen     = 20;
Wtildeinew = gen_Wtildei_12HEIplusenergy_child(recall,WtildeiS,beta,Xtildei1,Xtildei2,Utildei,n,iSigmae,Wistar,numgen);
Wtildei    = Wtildeinew;                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the MCMC traces and the posterior means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigmae_thin_trace = zeros(19,19,(nMCMC-nburn)/100);
Sigmau_thin_trace = zeros(19,19,(nMCMC-nburn)/100);
beta_thin_trace  = zeros(10,19,(nMCMC-nburn)/100);
Sigmae_mean_trace = zeros(19,19,(nMCMC-nburn)/100);
Sigmau_mean_trace = zeros(19,19,(nMCMC-nburn)/100);
beta_mean_trace  = zeros(10,19,(nMCMC-nburn)/100);

Sigmae_mean  = zeros(19,19);
Sigmau_mean  = zeros(19,19);
beta_mean    = zeros(10,19);

Wtildei_start = Wtildei;
kk=0;
sum_weightrecall = sum(weight.* recall);

disp('Start the MCMC');
for jjMCMC = 1:nMCMC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update W1 and W2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numgen     = 5;
Wtildeinew = gen_Wtildei_12HEIplusenergy_child(recall,Wtildei,beta,Xtildei1,Xtildei2,Utildei,n,iSigmae,Wistar,numgen);
Wtildei    = Wtildeinew;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate W-XB-U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt1 = zeros(n,19);
tt2 = zeros(n,19);
for jj = 1:19;
 tt1(:,jj) = (Xtildei1(:,:,jj) * beta(:,jj)) + Utildei(:,jj) ;
end;
for jj = 1:19;
 tt2(:,jj) = (Xtildei2(:,:,jj) * beta(:,jj)) + Utildei(:,jj) ;
end;
qqa = Wtildei(:,:,1) -  tt1;
qqb = (Wtildei(:,:,2) -  tt2) .*((recall==2)*ones(1,19));
qqaw = (weight*ones(1,19)).*(Wtildei(:,:,1) - tt1) ;
qqbw = ((recall==2).*weight*ones(1,19)).*(Wtildei(:,:,2) - tt2) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update iSigmae
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:5;
 rinew = updated_parameter_r_12HEIplusenergy_diffrecall(sum_weightrecall,i,r,theta,V,qqa,qqb,qqaw,qqbw,n);
 r(i,1)    =  rinew; 
 end;
for i = 1:25;
 thetainew = updated_parameter_theta_12HEIplusenergy(i,r,theta,V,qqa,qqb,qqaw,qqbw,n);
 theta(i,1)    =  thetainew; 
end;
for i = 2:2:12
  viinew = updated_parameter_Vii_12HEIplusenergy_diffrecall(sum_weightrecall,i,r,theta,V,qqa,qqb,qqaw,qqbw,n);
  V(i,i) = viinew;
end;
for i = 13:1:19
  viinew = updated_parameter_Vii_12HEIplusenergy_diffrecall(sum_weightrecall,i,r,theta,V,qqa,qqb,qqaw,qqbw,n);
  V(i,i) = viinew;
end;

for i = 4:2:12;
  for j = 1:(i-2);   
  vijnew = updated_parameter_Vij_12HEIplusenergy(i,j,r,theta,V,qqa,qqb,qqaw,qqbw,n);
  V(i,j) = vijnew;
  end;
end;
for i = 13:1:19;
  for j = 1:(i-1);   
  vijnew = updated_parameter_Vij_12HEIplusenergy(i,j,r,theta,V,qqa,qqb,qqaw,qqbw,n);
  V(i,j) = vijnew;
  end;
end;

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

Sigmae = V*V'; %'
iSigmae      = inv(Sigmae);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update iSigmaU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa         = ( (prior_Sigmau_doff-size(Sigmau,1)-1) .* prior_Sigmau_mean) + ( ((weight*ones(1,19)) .* Utildei)' * Utildei); %'
bb         = prior_Sigmau_doff + n;
aa  = aa + (.0011 .* (1 .* n) .* eye(19));
bb  = bb + (.0011 .* (1 .* n));
aa  = (aa + aa') ./ 2; %'

if min(eig(aa ./ bb)) <= 0.001;
    thecount2 = thecount2 + 1;
    if thecount2 == 1;
    disp(['Problem with Sigmau at step ',num2str(jjMCMC)]);
    disp('Sigmau');
    Sigmau
    disp('Sigmae');
    Sigmae
    disp('aa / bb');
    aa ./ bb;
    end;
end;
if min(eig(aa ./ bb)) > 0.001;
    Sigmau     = iwishrnd(aa,bb);
    iSigmau    = inv(Sigmau);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Utildei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qq1 = zeros(n,19);
qq2 = zeros(n,19);
for jj = 1:19;
     qq1(:,jj) = (Xtildei1(:,:,jj) * beta(:,jj)); 
end;
for jj = 1:19;
     qq2(:,jj) = (Xtildei2(:,:,jj) * beta(:,jj)); 
end;

ss1 = Wtildei(:,:,1) - qq1;
ss2 = (Wtildei(:,:,2) - qq2) .*((recall==2)*ones(1,19));
c1  = (ss1 + ss2) * iSigmae;
c2_2  = inv(iSigmau + (2.* iSigmae));
c2_1  = inv(iSigmau + iSigmae);
Utildei((1:n2),:) = (c1((1:n2),:)  * c2_2) + randn(n2,19) * sqrtm(c2_2);
Utildei((n2+1):n,:) = (c1((n2+1):n,:)  * c2_1) + randn(n1,19) * sqrtm(c2_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betanew = update_beta_12HEIplusenergy_child(recall,weight,iSigmau,iSigmae,prior_beta_cov, ...
                             Xtildei1,Xtildei2,Utildei,Wtildei,n,beta,prior_beta_mean);
beta    = betanew;

Sigmae_mean = Sigmae_mean + (Sigmae ./ nMCMC);
Sigmau_mean = Sigmau_mean + (Sigmau ./ nMCMC);
  beta_mean = beta_mean + (beta ./ nMCMC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get rid of the burn-in, Store every 100th results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if ((jjMCMC>nburn)& (rem(jjMCMC,100)==0));
 kk =kk+1; 
 Sigmae_mean_trace(:,:,kk) = Sigmae_mean;
 Sigmau_mean_trace(:,:,kk) = Sigmau_mean;
 for j= 1:19;
  beta_mean_trace(:,j,kk)  = beta_mean(:,j);
 end;
 Sigmae_thin_trace(:,:,kk) = Sigmae;
 Sigmau_thin_trace(:,:,kk) = Sigmau;
 for j= 1:19;
  beta_thin_trace(:,j,kk)  = beta(:,j);
 end;
 disp(num2str(jjMCMC));
 end;
end;
 
Sigmae_thin_trace = Sigmae_thin_trace .* (abs(Sigmae_thin_trace) > eps);

    save Test_Output/child_Sigmau_mean.mat Sigmau_mean_trace;
    save Test_Output/child_Sigmae_mean.mat Sigmae_mean_trace;
    save Test_Output/child_beta_mean.mat   beta_mean_trace;   

    save Test_Output/child_Sigmau.mat Sigmau_thin_trace;
    save Test_Output/child_Sigmae.mat Sigmae_thin_trace;
    save Test_Output/child_beta.mat   beta_thin_trace;   

 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Display the results
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   mm                = floor((nMCMC - nburn +eps) ./ 100);
%  xx = linspace(1,mm,mm)'; %'
%  figure(1);
%  k= 0 ;  
%  for i = 1:10;
%   for j = 1:10;
%      k = k+1;
%       subplot(10,10,k);
%       temp = zeros(mm,1);
%       for jj = 1:mm;
%        temp(jj,1) = Sigmae_thin_trace(i,j,jj);
%       end;  
%      plot(xx,temp);
%      title(['Sigmae(' num2str(i) ',' num2str(j) ')' ])
%    end;
%  end;
%  figure(2);
%  k= 0 ;  
%  for i = 1:10;
%   for j = 11:19;
%      k = k+1;
%       subplot(10,9,k);
%       temp = zeros(mm,1);
%       for jj = 1:mm;
%        temp(jj,1) = Sigmae_thin_trace(i,j,jj);
%       end;  
%      plot(xx,temp);
%      title(['Sigmae(' num2str(i) ',' num2str(j) ')' ])
%    end;
%  end;
%  figure(3);
%  k= 0 ;  
%  for i = 11:19;
%   for j = 11:19;
%      k = k+1;
%       subplot(9,9,k);
%       temp = zeros(mm,1);
%       for jj = 1:mm;
%        temp(jj,1) = Sigmae_thin_trace(i,j,jj);
%       end;  
%      plot(xx,temp);
%      title(['Sigmae(' num2str(i) ',' num2str(j) ')' ])
%    end;
%  end;
%  
%  figure(4);
%  k= 0 ;
%  for i = 1:10;
%   for j = 1:10;
%      k = k+1;
%       subplot(10,10,k);
%       temp = zeros(mm,1);
%       for jj = 1:mm;
%        temp(jj,1) = Sigmau_thin_trace(i,j,jj);
%       end;  
%      plot(xx,temp);
%      title(['Sigmau(' num2str(i) ',' num2str(j) ')' ])
%    end;
%  end;
%  figure(5);
%  k= 0 ;  
%  for i = 1:10;
%   for j = 11:19;
%      k = k+1;
%       subplot(10,9,k);
%       temp = zeros(mm,1);
%       for jj = 1:mm;
%        temp(jj,1) = Sigmau_thin_trace(i,j,jj);
%       end;  
%      plot(xx,temp);
%      title(['Sigmau(' num2str(i) ',' num2str(j) ')' ])
%    end;
%  end;
%  figure(6);
%  k= 0 ;  
%  for i = 11:19;
%   for j = 11:19;
%      k = k+1;
%       subplot(9,9,k);
%       temp = zeros(mm,1);
%       for jj = 1:mm;
%        temp(jj,1) = Sigmau_thin_trace(i,j,jj);
%       end;  
%      plot(xx,temp);
%      title(['Sigmau(' num2str(i) ',' num2str(j) ')' ])
%    end;
%  end;
%  
%      
%  figure(7);
%  k = 0 ;
%  for j = 1:10;
%    for i = 1:10;
%      k = k+1;
%      subplot(10,10,k);
%      temp = zeros(mm,1);
%      for jj = 1:mm;
%        temp(jj,1) = beta_thin_trace(i,j,jj);
%       end;  
%      plot(xx,temp);
%      title(['beta' num2str(j) '(' num2str(i) ')' ])
%    end;
%  end;
%  figure(8);
%  k = 0 ;
%  for j = 11:19;
%    for i = 1:10;
%      k = k+1;
%      subplot(9,10,k);
%      temp = zeros(mm,1);
%      for jj = 1:mm;
%        temp(jj,1) = beta_thin_trace(i,j,jj);
%       end;  
%      plot(xx,temp);
%      title(['beta' num2str(j) '(' num2str(i) ')' ])
%    end;
%  end;
  disp('*****************************************************************');
  disp(['number of simulation = ',num2str(mm)]);
  disp('*****************************************************************');
  disp(' ');
  format short;
  disp(' ');
  disp(' ');
  for i = 1:19;
    for j = i:19;
     disp(['$\Sigma_{\epsilon}(',num2str(i),',',num2str(j),')$&',num2str(mean(Sigmae_thin_trace(i,j,:))),'&',num2str(prctile(Sigmae_thin_trace(i,j,:),2.5)),'&',num2str(prctile(Sigmae_thin_trace(i,j,:),97.5)),'&',num2str(std(Sigmae_thin_trace(i,j,:))),'\\']);
    end;
  end;
  disp(' ');
  for i = 1:19;
    for j = i:19;
     disp(['$\Sigma_{u}(',num2str(i),',',num2str(j),')$&',num2str(mean(Sigmau_thin_trace(i,j,:))),'&',num2str(prctile(Sigmau_thin_trace(i,j,:),2.5)),'&',num2str(prctile(Sigmau_thin_trace(i,j,:),97.5)),'&',num2str(std(Sigmau_thin_trace(i,j,:))),'\\']);
    end;
  end;
  disp(' ');
  for j = 1:19;
    for i = 1:10;
     disp(['$\beta_{',num2str(j),'}(', num2str(i),',1)$&',num2str(mean(beta_thin_trace(i,j,:) )),'&',num2str(prctile(beta_thin_trace(i,j,:),2.5)),'&',num2str(prctile(beta_thin_trace(i,j,:),97.5)),'&',num2str(std(beta_thin_trace(i,j,:))),'\\']);
    end;
  end;
toc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After MCMC analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Test_Output/child_Sigmae.mat Sigmae_thin_trace;
load Test_Output/child_Sigmau.mat Sigmau_thin_trace;
load Test_Output/child_beta.mat beta_thin_trace;

a0 = zeros(13,1);
a0(1,1) = 0.5 .* min([HEI1_F_other(HEI1_F_other>0)' HEI2_F_other(HEI2_F_other>0)']');
a0(2,1) = 0.5 .* min([HEI1_F_WHL(HEI1_F_WHL>0)' HEI2_F_WHL(HEI2_F_WHL>0)']');

a0(3,1) = 0.5 .* min([HEI1_V_other(HEI1_V_other>0)' HEI2_V_other(HEI2_V_other>0)']');
a0(4,1) = 0.5 .* min([HEI1_V_DOL(HEI1_V_DOL>0)' HEI2_V_DOL(HEI2_V_DOL>0)']');

a0(5,1) = 0.5 .* min([HEI1_G_WHL(HEI1_G_WHL>0)' HEI2_G_WHL(HEI2_G_WHL>0)']');
a0(7,1) = 0.5 .* min([HEI1_G_other(HEI1_G_other>0)' HEI2_G_other(HEI2_G_other>0)']');

a0(6,1) = 0.5 .* min([HEI1_MILK(HEI1_MILK>0)'   HEI2_MILK(HEI2_MILK>0)']');

a0(8,1) = 0.5 .* min([HEI1_M_TOT(HEI1_M_TOT>0)' HEI2_M_TOT(HEI2_M_TOT>0)']');
a0(9,1) = 0.5 .* min([HEI1_OIL(HEI1_OIL>0)'     HEI2_OIL(HEI2_OIL>0)']');
a0(10,1)= 0.5 .* min([HEI1_SFAT(HEI1_SFAT>0)'   HEI2_SFAT(HEI2_SFAT>0)']');
a0(11,1)= 0.5 .* min([HEI1_SODI(HEI1_SODI>0)'   HEI2_SODI(HEI2_SODI>0)']');
a0(12,1)= 0.5 .* min([HEI1_SOFAAS_C(HEI1_SOFAAS_C>0)' HEI2_SOFAAS_C(HEI2_SOFAAS_C>0)']');
a0(13,1)= 0.5 .* min([energy1_other(energy1_other>0)' energy2_other(energy2_other>0)']'); %'


Sigmae_postmean = zeros(19,19);
Sigmau_postmean = zeros(19,19);
beta_postmean   = zeros(10,19);
for i = 1:19;
   for j = 1:19;
    Sigmae_postmean(i,j)= mean(Sigmae_thin_trace(i,j,:));
  end;
end;
for i = 1:19;
   for j = 1:19;
    Sigmau_postmean(i,j)= mean(Sigmau_thin_trace(i,j,:));
  end;
end;
for j = 1:19;
   for i = 1:10;
    beta_postmean(i,j)= mean(beta_thin_trace(i,j,:));
  end;
end;

B = 5;
T_bootstrap = zeros(n,13,B);
Xtildei          = zeros(n,10);

for b = 1:B;
U = randn(n,19) * sqrtm(Sigmau_postmean);
  
  for j =1:6; 
  Xtildei          = [ones(n,1) zeros(n,1) zeros(n,1) age sex race2 agesex agerace2 sexrace2 agesexrace2];
  intake           = (Xtildei(:,:) * beta_postmean(:,2*j) + U(:,2*j)); 
  min_intake       = sqrt(2)*(-1/lam(j,1)-mu(j,1))/sig(j,1);
  intake(intake < min_intake) = 0;
  temp             =  lam(j,1) .* (mu(j,1) + sig(j,1) ./ sqrt(2) .* intake) ;
  g_inv            = (1+temp) .^ (1/lam(j,1));
  der2_g_inv       = (sig(j,1)^2/2*(1-lam(j,1))).*((1+temp).^(-2+1/lam(j,1)));
  g_star           = g_inv + (0.5 * Sigmae_postmean(2*j,2*j)).* der2_g_inv;
  norm_cum         = normcdf(Xtildei(:,:) * beta_postmean(:,2*j-1) + U(:,2*j-1)); 
  T_wkday          = norm_cum .* g_star;
  T_wkday(T_wkday<a0(j,1)) = a0(j,1);

  Xtildei          = [ones(n,1) zeros(n,1) ones(n,1) age sex race2 agesex agerace2 sexrace2 agesexrace2];
  intake           = (Xtildei(:,:) * beta_postmean(:,2*j) + U(:,2*j)); 
  min_intake       = sqrt(2)*(-1/lam(j,1)-mu(j,1))/sig(j,1);
  intake(intake < min_intake) = 0;
  temp             =  lam(j,1) .* (mu(j,1) + sig(j,1) ./ sqrt(2) .* intake) ;
  g_inv            = (1+temp) .^ (1/lam(j,1));
  der2_g_inv       = (sig(j,1)^2/2*(1-lam(j,1))).*((1+temp).^(-2+1/lam(j,1)));
  g_star           = g_inv + (0.5 * Sigmae_postmean(2*j,2*j)).* der2_g_inv;
  norm_cum         = normcdf(Xtildei(:,:) * beta_postmean(:,2*j-1) + U(:,2*j-1)); 
  T_wkend              = norm_cum .* g_star;
  T_wkend(T_wkend<a0(j,1)) = a0(j,1);
  
  T_bootstrap(:,j,b) = (4/7)*T_wkday + (3/7)*T_wkend;
  end;

  for j= 7:13;
  Xtildei          = [ones(n,1) zeros(n,1) zeros(n,1) age sex race2 agesex agerace2 sexrace2 agesexrace2];
  intake           = (Xtildei(:,:) * beta_postmean(:,j+6) + U(:,j+6)); 
  min_intake       = sqrt(2)*(-1/lam(j,1)-mu(j,1))/sig(j,1);
  intake(intake < min_intake) = 0;
  temp             =  lam(j,1) .* (mu(j,1)+ sig(j,1) ./ sqrt(2) .* intake) ;
  g_inv            = (1+temp) .^ (1/lam(j,1));
  der2_g_inv       = (((sig(j,1)^2)/2)*(1-lam(j,1))).*((1+temp).^(-2 + 1/lam(j,1)));
  g_star           = g_inv + (0.5 * Sigmae_postmean(j+6,j+6)).* der2_g_inv;
  T_wkday              = g_star;
  T_wkday(T_wkday<a0(j,1)) = a0(j,1);
    
  Xtildei          = [ones(n,1) zeros(n,1) ones(n,1) age sex race2 agesex agerace2 sexrace2 agesexrace2];
  intake           = (Xtildei(:,:) * beta_postmean(:,j+6) + U(:,j+6)); 
  min_intake       = sqrt(2)*(-1/lam(j,1)-mu(j,1))/sig(j,1);
  intake(intake < min_intake) = 0;
  temp             =  lam(j,1) .* (mu(j,1)+ sig(j,1) ./ sqrt(2) .* intake) ;
  g_inv            = (1+temp) .^ (1/lam(j,1));
  der2_g_inv       = (((sig(j,1)^2)/2)*(1-lam(j,1))).*((1+temp).^(-2 + 1/lam(j,1)));
  g_star           = g_inv + (0.5 * Sigmae_postmean(j+6,j+6)).* der2_g_inv;
  T_wkend              = g_star;
  T_wkend(T_wkend <a0(j,1)) = a0(j,1);

  T_bootstrap(:,j,b)= (4/7)*T_wkday + (3/7)*T_wkend;
  end;

end;

% order:  F_other, F_WHL, V_other, V_DOL, G_WHL, MILK, G_other, M_TOT, OIL, SFAT, SODI, SOFAAS, energy_other

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%  USUAL INTAKE
 q1  = zeros(101,13);
 TT_reorder = zeros(n,13,B);
 TT_reorder= T_bootstrap;
 TT_reorder(:,1,:) = T_bootstrap(:,1,:)+T_bootstrap(:,2,:);
 TT_reorder(:,5,:) = T_bootstrap(:,3,:)+T_bootstrap(:,4,:);
 TT_reorder(:,3,:) = T_bootstrap(:,5,:)+T_bootstrap(:,7,:);
 TT_reorder(:,4,:) = T_bootstrap(:,5,:);
 TT_reorder(:,6,:) = T_bootstrap(:,4,:);
 TT_reorder(:,7,:) = T_bootstrap(:,6,:);
 TT_reorder(:,13,:)= T_bootstrap(:,13,:) +9.* T_bootstrap(:,10,:)+T_bootstrap(:,12,:);

% order:  F_TOT, F_WHL, G_TOT, G_WHL, V_TOT, V_DOL, MILK,  M_TOT, OIL, SFAT, SODI, SOFAAS, energy

% R= zeros(13,13);
%  for i =1:13;
%  for j =i:13; 
%    S1 = reshape(TT_reorder(:,i,:),n*B,1);
%    S2 = reshape(TT_reorder(:,j,:),n*B,1);
%    R(i,j) = corr(S1,S2);
%  end;
%  end;
%  R=R+R'-eye(13); %'
% format bank;
%  R = floor(100 .* (R+0.001)) ./ 100
%%  save usualintake_R.txt R -ascii;
% for  i= 1:13;
% disp(['&',num2str(R(i,1)),'&',num2str(R(i,2)),'&',num2str(R(i,3)),'&',num2str(R(i,4)),'&',num2str(R(i,5)),'&',num2str(R(i,6)),'&',num2str(R(i,7)),'&',num2str(R(i,8)),'&',num2str(R(i,9)),'&',num2str(R(i,10)),'&',num2str(R(i,11)),'&',num2str(R(i,12)),'&',num2str(R(i,13)),'\\']);
% end;

  for j =1:13;
  T = reshape(TT_reorder(:,j,:),n,B);
  x = linspace(0, max(max(T)),max(200,10*max(max(T))));
  K = size(x,2);
  tilde = x(2)-x(1);
   %mean
  q1(100,j) = sum( weight .* sum(T')')/(B*n);
   %sd
  q1(101,j) = sqrt(sum( weight .* sum(T.^2')')/(B*n)-q1(100,j)^2);
   c = zeros(K,1);
  for i=1:K;
  c(i) = sum(weight.*(sum((T<=x(i))')'))/(B*n);
  end;
  p = linspace(0.01,0.99,99);
  for i=1:99; 
  index=max(find(c<p(i)));
  q1(i,j)=  (p(i)-c(index))*(x(index+1)-x(index))/(c(index+1)-c(index)) + x(index);
  end;
  end;
  save Test_Output/q1.mat q1;
  for j = 1:13;
  disp(['&',num2str(q1(100,j)),'&',num2str(q1(101,j)),'&', num2str(q1(5,j)),'&',num2str(q1(10,j)),'&',num2str(q1(25,j)),'&',num2str(q1(50,j)),'&',num2str(q1(75,j)),'&',num2str(q1(90,j)),'&',num2str(q1(95,j)),'\\']);
  end;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%  adjusted INTAKE
  q2  = zeros(101,12);
 energy = TT_reorder(:,13,:);
   T_ad = zeros(n,12,B);
   for i =1:9;
    T_ad(:,i,:) = 1000.* reshape(TT_reorder(:,i,:)./energy,n,B);
   end;
   i=10;
   T_ad(:,10,:) = 900 .* reshape(TT_reorder(:,i,:)./energy,n,B);
   i=11;
    T_ad(:,11,:) = reshape(TT_reorder(:,i,:)./energy,n,B);
   i=12;
    T_ad(:,12,:) = 100.* reshape(TT_reorder(:,i,:)./energy,n,B); 

    R= zeros(12,12);
    for i =1:12;
    for j =i:12; 
      S1 = reshape(T_ad(:,i,:),n*B,1);
      S2 = reshape(T_ad(:,j,:),n*B,1);
      R(i,j) = corr(S1,S2);
    end;
    end;
   R=R+R'-eye(12); %'
   format bank;
   R = floor(100 .* (R+0.001)) ./ 100
   save Test_Output/adjustedintake_R.txt R -ascii;
   for  i= 1:12;
   disp(['&',num2str(R(i,1)),'&',num2str(R(i,2)),'&',num2str(R(i,3)),'&',num2str(R(i,4)),'&',num2str(R(i,5)),'&',num2str(R(i,6)),'&',num2str(R(i,7)),'&',num2str(R(i,8)),'&',num2str(R(i,9)),'&',num2str(R(i,10)),'&',num2str(R(i,11)),'&',num2str(R(i,12)),'\\']);
   end;
 
 for j =1:12;
 T = reshape(T_ad(:,j,:),n,B);
 x = linspace(0, max(max(T)),max(200,10*max(max(T))));
 K = size(x,2);
 tilde = x(2)-x(1);
 %mean
 q2(100,j)= sum( weight .* sum(T')')/(B*n);
 %sd
 q2(101,j) = sqrt(sum( weight .* sum(T.^2')')/(B*n)-q2(100,j)^2);
 c = zeros(K,1);
 for i=1:K;
 c(i) = sum(weight.*(sum((T<=x(i))')'))/(B*n);
 end;
 p = linspace(0.01,0.99,99);
 for i=1:99; 
 index=max(find(c<p(i)));
 q2(i,j)=  (p(i)-c(index))*(x(index+1)-x(index))/(c(index+1)-c(index)) + x(index);
 end;
end;
 save Test_Output/q2.mat q2;
 for j = 1:12;
 disp(['&',num2str(q2(100,j)),'&',num2str(q2(101,j)),'&', num2str(q2(5,j)),'&',num2str(q2(10,j)),'&',num2str(q2(25,j)),'&',num2str(q2(50,j)),'&',num2str(q2(75,j)),'&',num2str(q2(90,j)), '&',num2str(q2(95,j)),'&',num2str(q2(99,j)),'\\']);
 end;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%  score
 %order:  F_TOT, F_WHL, G_TOT, G_WHL,V_TOT, V_DOL, MILK, M_TOT, OIL, SFAT, SODI, SOFAAS, Total score
 
 q3 = zeros(101,13);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The realization of the variables in the HEI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 KCAL  = T_bootstrap(:,13,:) +9.* T_bootstrap(:,10,:)+T_bootstrap(:,12,:); % energy;
 F_TOT = T_bootstrap(:,1,:)+T_bootstrap(:,2,:); 
 F_WHL = T_bootstrap(:,2,:); 
 V_TOT = T_bootstrap(:,3,:)+T_bootstrap(:,4,:); 
 V_DOL = T_bootstrap(:,4,:); 
 G_WHL = T_bootstrap(:,5,:);
 MILK  = T_bootstrap(:,6,:);  
 G_TOT = T_bootstrap(:,5,:)+T_bootstrap(:,7,:);  
 M_TOT = T_bootstrap(:,8,:);
 OIL   = T_bootstrap(:,9,:);
 SFAT  = T_bootstrap(:,10,:);
 SODI  = T_bootstrap(:,11,:);
 SOFAAS_C =  T_bootstrap(:,12,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form the quantiles of usual intakes. The first 99 are the
% quantiles, the 100th is the mean, and the 101st is the standard deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:13;
 S = reshape(T_bootstrap(:,j,:),n,B);
 x = linspace(0, max(max(S)),max(200,10*max(max(S))));
 K = size(x,2);
 tilde = x(2)-x(1);
 %mean
 q3(100,j) = sum( weight .* sum(S')')/(B*n);
 %sd
 q3(101,j) = sqrt(sum( weight .* sum(S.^2')')/(B*n)-q3(100,j)^2);
 c = zeros(K,1);
 for i=1:K;
 c(i) = sum(weight.*(sum((S<=x(i))')'))/(B*n);
 end;
 p = linspace(0.01,0.99,99);
 for i=1:99; 
  index=max(find(c<p(i)));
  if size(index,1)==0;
   q3(i,j)= 0;
 end;
  if size(index,1)~=0;
   q3(i,j) =  (p(i)-c(index))*(x(index+1)-x(index))/(c(index+1)-c(index)) + x(index);
  end;
 end;
end;
save Test_Output/quantiles_Usual_Intakes.mat q3;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The realization of the energy adjusted densities in the HEI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Density_bootstrap = zeros(n,12,B);
 DEN1 =  F_TOT./(KCAL./1000);
 DEN2 =  F_WHL./(KCAL./1000);
 DEN3 =  G_TOT./(KCAL./1000);
 DEN4 =  G_WHL./(KCAL./1000);
 DEN5 =  V_TOT./(KCAL./1000);
 DEN6 =  V_DOL./(KCAL./1000);
 DEN7 =  MILK./(KCAL./1000);
 DEN8 =  M_TOT./(KCAL./1000);
 DEN9 =  OIL./(KCAL./1000);
 DEN10 = 100*(SFAT*9)./KCAL;
 DEN11 = SODI./(KCAL./1000);
 DEN12 = 100.*(SOFAAS_C./KCAL);
Density_bootstrap = zeros(n,12,B);
Density_bootstrap(:,1,:) = DEN1; 
Density_bootstrap(:,2,:) = DEN2; 
Density_bootstrap(:,3,:) = DEN3; 
Density_bootstrap(:,4,:) = DEN4; 
Density_bootstrap(:,5,:) = DEN5; 
Density_bootstrap(:,6,:) = DEN6; 
Density_bootstrap(:,7,:) = DEN7; 
Density_bootstrap(:,8,:) = DEN8; 
Density_bootstrap(:,9,:) = DEN9; 
Density_bootstrap(:,10,:) = DEN10; 
Density_bootstrap(:,11,:) = DEN11; 
Density_bootstrap(:,12,:) = DEN12; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form the quantiles of usual intakes. The first 99 are the
% quantiles, the 100th is the mean, and the 101st is the standard deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:12;
 S = reshape(Density_bootstrap(:,j,:),n,B);
 x = linspace(0, max(max(S)),max(200,10*max(max(S))));
 K = size(x,2);
 tilde = x(2)-x(1);
 %mean
 q3(100,j) = sum( weight .* sum(S')')/(B*n);
 %sd
 q3(101,j) = sqrt(sum( weight .* sum(S.^2')')/(B*n)-q3(100,j)^2);
 c = zeros(K,1);
 for i=1:K;
 c(i) = sum(weight.*(sum((S<=x(i))')'))/(B*n);
 end;
 p = linspace(0.01,0.99,99);
 for i=1:99; 
  index=max(find(c<p(i)));
  if size(index,1)==0;
   q3(i,j)= 0;
 end;
  if size(index,1)~=0;
   q3(i,j) =  (p(i)-c(index))*(x(index+1)-x(index))/(c(index+1)-c(index)) + x(index);
  end;
 end;
end;
save Test_Output/quantiles_Usual_Densities.mat q3;

 
 
 
 HEI =zeros(n,13,B);
 HEI(:,1,:) = min(5,5.*(DEN1./0.8));
 HEI(:,2,:) = min(5,5.*(DEN2./0.4));
 HEI(:,3,:) =min(5, 5.*(DEN3./3));
 HEI(:,4,:) =min(5, 5.*(DEN4./1.5));
 HEI(:,5,:) =min(5, 5.*(DEN5./1.1));
 HEI(:,6,:) =min(5, 5.*(DEN6./0.4));
 HEI(:,7,:) =min(10, 10.*(DEN7./1.3));
 HEI(:,8,:) =min(10, 10.*(DEN8./2.5));
 HEI(:,9,:) =min(10, 10.*(DEN9./12));
 HEI(:,10,:) = 0.*(DEN10>=15)  + 10.*(DEN10<=7)  +(8 - ( 8 .* (DEN10-10)./5 )).*((DEN10>10)&(DEN10<15))+(10 - (2 .* (DEN10-7)./3 )).*((DEN10>7)&(DEN10<=10));
 HEI(:,11,:) = 0.*(DEN11>=2000)+ 10.*(DEN11<=700)+(8 - ( 8 .* (DEN11-1100)./(2000-1100) )).*((DEN11>=1100)&(DEN11<2000))+(10 - (2 .* (DEN11-700)./(1100-700) )).*((DEN11>700)&(DEN11<1100));
 HEI(:,12,:) = 0.*(DEN12>=50)  + 20.*(DEN12<=20) +(20-(20 .* (DEN12-20)./(50 - 20))).*((DEN12>20)&(DEN12<50));
 HEI(:,13,:) = HEI(:,1,:)+HEI(:,2,:)+HEI(:,3,:)+HEI(:,4,:)+HEI(:,5,:)+HEI(:,6,:)+HEI(:,7,:)+HEI(:,8,:)+HEI(:,9,:)+HEI(:,10,:)+HEI(:,11,:)+HEI(:,12,:);    

  R= zeros(12,1);
  for i =1:12;
    S1 = reshape(HEI(:,i,:),n*B,1);
    S2 = reshape(HEI(:,13,:),n*B,1)-reshape(HEI(:,i,:),n*B,1);
    R(i,1) = corr(S1,S2);
  end;
 R
 save Test_Output/score_R_model.mat R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form the quantiles of usual HEI Total Score. The first 99 are the
% quantiles, the 100th is the mean, and the 101st is the standard deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:13;
 S = reshape(HEI(:,j,:),n,B);
 x = linspace(0, max(max(S)),max(200,10*max(max(S))));
 K = size(x,2);
 tilde = x(2)-x(1);
 %mean
 q3(100,j) = sum( weight .* sum(S')')/(B*n);
 %sd
 q3(101,j) = sqrt(sum( weight .* sum(S.^2')')/(B*n)-q3(100,j)^2);
 c = zeros(K,1);
 for i=1:K;
 c(i) = sum(weight.*(sum((S<=x(i))')'))/(B*n);
 end;
 p = linspace(0.01,0.99,99);
 for i=1:99; 
  index=max(find(c<p(i)));
  if size(index,1)==0;
   q3(i,j)= 0;
 end;
  if size(index,1)~=0;
   q3(i,j) =  (p(i)-c(index))*(x(index+1)-x(index))/(c(index+1)-c(index)) + x(index);
  end;
 end;
end;

save Test_Output/quantiles_Usual_HEI_Scores.mat q3;

save Test_Output/q3.mat q3;
 for j = 1:13;
 disp(['&',num2str(q3(100,j)),'&',num2str(q3(101,j)),'&', num2str(q3(5,j)),'&',num2str(q3(10,j)),'&',num2str(q3(25,j)),'&',num2str(q3(50,j)),'&',num2str(q3(75,j)),'&',num2str(q3(90,j)) ,'&',num2str(q3(95,j)),'&',num2str(q3(99,j)),'\\']);
 end;

% Figure 1 in paper: the estimated percentiles of the total score
figure(9);
load Test_Output/q3.mat q3;
h2= plot(p,q3(1:99,13), 'b-');
set(h2,'LineWidth',3);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dietary consistency 
 load Test_Output/q3.mat q3;
 k=[q3(50,1) q3(50,2) q3(50,5) q3(50,6) q3(50,4) q3(50,7)];
 HEI_epi = [HEI(:,1,:), HEI(:,2,:), HEI(:,5,:), HEI(:,6,:),HEI(:,4,:), HEI(:,7,:)];
 HEI_epi1 = reshape(HEI_epi(:,1,:),n,B)>=k(1) ;
 HEI_epi2 = reshape(HEI_epi(:,2,:),n,B)>=k(2) ;
 HEI_epi3 = reshape(HEI_epi(:,3,:),n,B)>=k(3) ;
 HEI_epi4 = reshape(HEI_epi(:,4,:),n,B)>=k(4) ;
 HEI_epi5 = reshape(HEI_epi(:,5,:),n,B)>=k(5) ;
 HEI_epi6 = reshape(HEI_epi(:,6,:),n,B)>=k(6) ;
 epi_product = HEI_epi1.* HEI_epi2.* HEI_epi3.* HEI_epi4 .* HEI_epi5.* HEI_epi6;
 cdf = sum(sum((weight)* ones(1,B).*epi_product))/(n*B);
 totalp =0;
 for i = 1:n;
  for b= 1:B;
  totalp = totalp + weight(i) .* HEI_epi1(i,b) .*HEI_epi2(i,b).*HEI_epi3(i,b) .*HEI_epi4(i,b) .*HEI_epi5(i,b) .*HEI_epi6(i,b); 
 end;
 end;
 disp(['The percentage of child aged 2 to 8 with HEI score exceeds the median HEI score on all 6 episodically consumed foods:',num2str(100.* totalp./(n.*B)),'%']);

 load Test_Output/q3.mat q3;
 k=q3(50,1:12);
 HEI_epi1 = reshape(HEI(:,1,:),n,B)>=k(1) ;
 HEI_epi2 = reshape(HEI(:,2,:),n,B)>=k(2) ;
 HEI_epi3 = reshape(HEI(:,3,:),n,B)>=k(3) ;
 HEI_epi4 = reshape(HEI(:,4,:),n,B)>=k(4) ;
 HEI_epi5 = reshape(HEI(:,5,:),n,B)>=k(5) ;
 HEI_epi6 = reshape(HEI(:,6,:),n,B)>=k(6) ;
 HEI_epi7 = reshape(HEI(:,7,:),n,B)>=k(7) ;
 HEI_epi8 = reshape(HEI(:,8,:),n,B)>=k(8) ;
 HEI_epi9 = reshape(HEI(:,9,:),n,B)>=k(9) ;
 HEI_epi10 = reshape(HEI(:,10,:),n,B)>=k(10) ;
 HEI_epi11 = reshape(HEI(:,11,:),n,B)>=k(11) ;
 HEI_epi12 = reshape(HEI(:,12,:),n,B)>=k(12) ;
 epi_product = HEI_epi1.* HEI_epi2.* HEI_epi3.* HEI_epi4 .* HEI_epi5.* HEI_epi6.* HEI_epi7.* HEI_epi8.* HEI_epi9.* HEI_epi10.* HEI_epi11.* HEI_epi12;
 cdf = sum(sum((weight)* ones(1,B).*epi_product))/(n*B);
 disp(['The percentage of child aged 2 to 8 with HEI score exceeds the median HEI score on all 12 HEI components:',num2str(100.* cdf),'%']);

 load Test_Output/q3.mat q3;
 k=q3(1:99,1:12);
 cdf = zeros(99,1);
 for i = 1: size(k,1);
 HEI_epi1 = reshape(HEI(:,1,:),n,B)>=k(i,1) ;
 HEI_epi2 = reshape(HEI(:,2,:),n,B)>=k(i,2) ;
 HEI_epi3 = reshape(HEI(:,3,:),n,B)>=k(i,3) ;
 HEI_epi4 = reshape(HEI(:,4,:),n,B)>=k(i,4) ;
 HEI_epi5 = reshape(HEI(:,5,:),n,B)>=k(i,5) ;
 HEI_epi6 = reshape(HEI(:,6,:),n,B)>=k(i,6) ;
 HEI_epi7 = reshape(HEI(:,7,:),n,B)>=k(i,7) ;
 HEI_epi8 = reshape(HEI(:,8,:),n,B)>=k(i,8) ;
 HEI_epi9 = reshape(HEI(:,9,:),n,B)>=k(i,9) ;
 HEI_epi10 = reshape(HEI(:,10,:),n,B)>=k(i,10) ;
 HEI_epi11 = reshape(HEI(:,11,:),n,B)>=k(i,11) ;
 HEI_epi12 = reshape(HEI(:,12,:),n,B)>=k(i,12) ;
 epi_product = HEI_epi1.* HEI_epi2.* HEI_epi3.* HEI_epi4 .* HEI_epi5.* HEI_epi6.* HEI_epi7.* HEI_epi8.* HEI_epi9.* HEI_epi10.* HEI_epi11.* HEI_epi12;
 cdf(i) = sum(sum((weight)* ones(1,B).*epi_product))/(n*B);
 end;

%Figure 3 in paper
figure(10);
h2= plot(p.*100, cdf, 'b-');
set(h2,'LineWidth',3);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% conditional 
q4 = zeros(101,12);
totalscorelt50 =(reshape(HEI(:,13,:),n,B) <=50); 
denominator = (sum(weight .* sum(totalscorelt50')'));   

for j =[2,4,6,12];
 T = reshape(T_ad(:,j,:),n,B);
 x = linspace(0, max(max(T)),max(200,10*max(max(T))));
 K = size(x,2);
 tilde = x(2)-x(1);
 %mean
 q4(100,j)= sum( weight .* sum((T.* totalscorelt50)')') / denominator;
 %sd
 q4(101,j) = sqrt(sum( weight .* sum((T.* totalscorelt50).^2')')/denominator-q4(100,j)^2);
 
 c = zeros(K,1);
 for i=1:K;
 c(i) = sum(weight.*(sum(((T<=x(i)).*totalscorelt50)')'))/denominator;
 end;
 p = linspace(0.01,0.99,99);
 for i=1:99; 
  index=max(find(c<p(i)));
  if size(index,1)==0;
    q4(i,j)= 0;
  end;
  if size(index,1)~=0;
  q4(i,j)=  (p(i)-c(index))*(x(index+1)-x(index))/(c(index+1)-c(index)) + x(index);
  end;
 end;
end;
 save Test_Output/q4_lt50.mat q4;
  for j =[2,4,6,12];
  disp(['&',num2str(q4(100,j)),'&',num2str(q4(101,j)),'&', num2str(q4(5,j)),'&',num2str(q4(10,j)),'&',num2str(q4(25,j)),'&',num2str(q4(50,j)),'&',num2str(q4(75,j)),'&',num2str(q4(90,j)), '&',num2str(q4(95,j)),'\\']);
  end;

q4 = zeros(101,12);
totalscoregt50 =(reshape(HEI(:,13,:),n,B) >50); 
denominator = (sum(weight .* sum(totalscoregt50')'));

for j =[2,4,6,12];
 T = reshape(T_ad(:,j,:),n,B);
 x = linspace(0, max(max(T)),max(200,10*max(max(T))));
 K = size(x,2);
 tilde = x(2)-x(1);
 %mean
 q4(100,j)= sum( weight .* sum((T.* totalscoregt50)')') / denominator;
 %sd
 q4(101,j) = sqrt(sum( weight .* sum((T.* totalscoregt50).^2')')/denominator-q4(100,j)^2);
 
 c = zeros(K,1);
 for i=1:K;
 c(i) = sum(weight.*(sum(((T<=x(i)).*totalscoregt50)')'))/denominator;
 end;
 p = linspace(0.01,0.99,99);
 for i=1:99; 
  index=max(find(c<p(i)));
  if size(index,1)==0;
    q4(i,j)= 0;
  end;
  if size(index,1)~=0;
  q4(i,j)=  (p(i)-c(index))*(x(index+1)-x(index))/(c(index+1)-c(index)) + x(index);
  end;
 end;
end;
save Test_Output/q4_gt50.mat q4;
for j =[2,4,6,12];
disp(['&',num2str(q4(100,j)),'&',num2str(q4(101,j)),'&', num2str(q4(5,j)),'&',num2str(q4(10,j)),'&',num2str(q4(25,j)),'&',num2str(q4(50,j)),'&',num2str(q4(75,j)),'&',num2str(q4(90,j)), '&',num2str(q4(95,j)),'\\']);
end;

%Figure 2 in paper
load Test_Output/q4_lt50.mat q4;
q4lt50=q4;
load Test_Output/q4_gt50.mat q4;
q4gt50=q4;
p = linspace(0.01,0.99,99);

 subplot(2,2,1)
 h2= plot(p, q4lt50(1:99,2), 'b-', p,q4gt50(1:99,2), 'r--');
 set(h2,'LineWidth',3);
 S=['Whole Fruit']
 h=title(S);
 set(h,'FontWeight','bold')
 set(h,'FontSize',14)

 subplot(2,2,2)
 h2= plot(p, q4lt50(1:99,6), 'b-', p,q4gt50(1:99,6), 'r--');
 set(h2,'LineWidth',3);
 S=['Whole Grains']
 h=title(S);
 set(h,'FontWeight','bold')
 set(h,'FontSize',14)

subplot(2,2,3)
 h2= plot(p, q4lt50(1:99,4), 'b-', p,q4gt50(1:99,4), 'r--');
 set(h2,'LineWidth',3);
 S=['DOL']
 h=title(S);
 set(h,'FontWeight','bold')
 set(h,'FontSize',14)

 subplot(2,2,4)
 h2= plot(p, q4lt50(1:99,12), 'b-', p,q4gt50(1:99,12), 'r--');
 set(h2,'LineWidth',3);
 S=['SoFAAS']
 h=title(S);
 set(h,'FontWeight','bold')
 set(h,'FontSize',14)

diary off;


