%Psychometric function with logistic regression
%Input: X,XMEAN,alpha,beta
%Output: estimated response data (proportion of one of two alternative
%responses)
function [rprop] = psychfunc(X,XMEAN,alpha,beta)
rprop=zeros(1,length(X));
for i=1:length(X)
    rprop(i) = exp(alpha + beta*(X(i) - XMEAN))/(1+exp(alpha + beta*(X(i) - XMEAN)));
end
end