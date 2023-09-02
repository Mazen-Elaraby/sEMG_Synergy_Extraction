function [KL] = KL_Divergence(org,rec)

NoC = size(org,2);

KL_arr = zeros(1,NoC);

%edge cases:
    %0/0 divisions -> NaN
    %0/x division  -> log2 returns -Inf
    %x/0 division  -> Inf

%adding softening coefficent
s_k = 0.00001; %all signal values are in the order of 4 decimals 

for i=1:NoC
    %Probability Density of Target distribution
    p = histcounts(org(:,i),'Normalization','probability');
    %Probability of Approximate Distribution
    q = histcounts(rec(:,i),length(p),'Normalization','probability');
    
    q = q + s_k; 
    
    KL_arr(i) = sum(p .* log2( (p./q) + s_k) );
end

KL = mean(KL_arr);

end