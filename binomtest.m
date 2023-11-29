function pv = binomtest(pvs,alpha)
k = length(pvs);
% compute test statistic
statistic = sum(pvs <= alpha);
pv = sum(binopdf(statistic:k,k,alpha));
end