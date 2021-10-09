function genww = gen_truncated_normals(trunc_value,startxi,numgen);
%
% This generates standard normals truncated from the left
% at trunc_value, with starting values startxi,
% using a rejection sampler devised by C. P. Robert
% Statistics and Computing, 1995, volumn 5, pp 121-125.
%
% The rejection sampler is only used for truncation values
% > 0, because it is remarkably bad for truncation values
% < 0, i.e., it tends to stay where it starts in this case.
%
% INPUT:
%      trunc_value: the random variable is truncated
%                   at the left from trunc_value, a vector
%      startxi:     starting values
%      numgen:      number of times you try (recommended = 50)
%
% OUTPUT:
%      truncated normals, the same dimension as trunc_value
%
%
n      = size(trunc_value,1);
alpha  = (trunc_value + sqrt(4 + (trunc_value .^ 2))) ./ 2;
thesign = (trunc_value >= 0); % Whether the truncation point is positive
genww  = trunc_value .* (trunc_value > 0);
temp2  = randn(n,1);
for jj = 1:numgen;
xicand = trunc_value - ( (1 ./ alpha) .* log(rand(n,1)));
mmmm   = (rand(n,1) < exp(-.5 .* ( (xicand - alpha) .^ 2)));
temp1  = (xicand  .* (mmmm == 1)) + (genww .* (mmmm == 0));
ssss   = randn(n,1);
temp2  = (temp2 .* (ssss < trunc_value)) + ...
               (ssss .* (ssss >= trunc_value));
genww  = (temp2 .* (thesign == 0)) + (temp1 .* (thesign == 1));
end;
genww  = (genww .* (genww > trunc_value)) + ((trunc_value + eps) .* (genww <= trunc_value));

  

