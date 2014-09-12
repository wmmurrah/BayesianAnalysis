%Inverse psychometric function with logistic regression
%Input: Y,alpha,beta
%Output: estimated stimulus intensity
function [stim] = psychfunc_inv(Y,XMEAN,alpha,beta)
stim=zeros(1,length(alpha));
for i=1:length(alpha)
    stim(i) = (log(-Y/(Y-1))-alpha(i))/beta(i) + XMEAN;
end
end