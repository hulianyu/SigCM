function pval_samples = SigCM_intra(X,target_cluster)
Alpha = 0.05;
Dk = X(target_cluster,:);
pvals= sample_Fisher_exact_test(X, Dk);
numClu = length(target_cluster);
pval_samples = zeros(numClu,1);
parfor i=1:numClu
    pval_samples(i,1) = binomtest(pvals(i,:),Alpha);
end
pval_samples = [target_cluster pval_samples];
end

function p = sample_Fisher_exact_test(X, X_target_cluster)
% He Z, Zhao C, Liang H, et al. Protein complexes identification with family-wise error rate control[J].
% IEEE/ACM transactions on computational biology and bioinformatics, 2019, 17(6): 2062-2073.
% Lancichinetti A, Radicchi F, Ramasco J J, et al. Finding statistically significant communities in networks[J]. 
% PloS one, 2011, 6(4): e18961.
% record of revisions:
%     date               programmer              description of change
% -----------        -----------------          ------------------------
%  Oct 11, 2023           Lianyu Hu              Original code version
%% For each attribute of a cluster, compute p-value on point s
[N,num_attr] = size(X);
Nk = size(X_target_cluster,1);
p = zeros(Nk,num_attr);
parfor attr=1:num_attr
    X_attr = X(:,attr);
    [value, tab] = histRate(X_attr);
    for th = 1:Nk
        s = X_target_cluster(th,attr);
        k_i = tab(value==s);
        k_in =  sum(X_target_cluster(:,attr)==s);
        k_out = k_i - k_in;
        volume = min(Nk-k_in, k_out);
        p(th,attr) = sum(hygepdf(k_in:k_in+volume,N,Nk,k_i));
    end
end
end

function pval = binomtest(pvals,Alpha)
M = length(pvals);
% compute test statistic
statistic = sum(pvals <= Alpha);
pval = sum(binopdf(statistic:M,M,Alpha));
end