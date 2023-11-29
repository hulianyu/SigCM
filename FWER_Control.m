function [signum,sigp,siglist] = FWER_Control(pvs)
%% FWER: Bonferroni method
alpha = 0.05;
m = size(pvs,1);
adjusted_alpha = alpha / m; % adjusted significance level
siglist = find(pvs(:,2) <= adjusted_alpha);
if isempty(siglist)
    siglist = [];
    signum = 0;
    sigp = 0;
else
    siglist = pvs(siglist,:);
    signum = length(siglist);
    sigp = signum/m;
end
end