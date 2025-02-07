function pv_point = Assign_Control(sList,X,k,pi)
% Calculate the percentage of samples controlled by the assignment
%% assignment: pval(thisK)<pval(otherK)
OtherClusters = setdiff(unique(pi),k);
numClu = size(sList,1);
numOther = length(OtherClusters);
pv_point = zeros(numClu,numOther);
for kth = 1:length(OtherClusters)
otherK = pi==OtherClusters(kth);
OtherS = X(otherK,:);
% X_s = X(otherK,:);
pvs = Other_single_point_fisher_exactG(X, sList, OtherS);
for sth=1:numClu
    pv_point(sth,kth) = binomtest(pvs(sth,:),0.05);
end
end
end

function p = Other_single_point_fisher_exactG(X, sList, OtherS)
% He Z, Zhao C, Liang H, et al. Protein complexes identification with family-wise error rate control[J].
% IEEE/ACM transactions on computational biology and bioinformatics, 2019, 17(6): 2062-2073.
% Lancichinetti A, Radicchi F, Ramasco J J, et al. Finding statistically significant communities in networks[J]. 
% PloS one, 2011, 6(4): e18961.
% record of revisions:
%     date               programmer              description of change
% -----------        -----------------          ------------------------
%  Nov 20, 2023           Lianyu Hu              Original code version
%% for each attribute value u of other cluster OtherS, compute p-value on point sList
[N,num_attr] = size(X);
i_S = size(sList,1);
i_X = X(sList(:,1),:);
D_S = size(OtherS,1);
p = zeros(i_S,num_attr);
for attr=1:num_attr
    X_attr = X(:,attr);
    [value, tab] = histRate(X_attr);
    for th = 1:i_S
        s = i_X(th,attr);
        k_i = tab(value==s);
        k_in =  sum(OtherS(:,attr)==s);
        k_out = k_i - k_in;
        volume = min(D_S-k_in, k_out);
        p(th,attr) = sum(hygepdf(k_in:k_in+volume,N,D_S,k_i));
    end
end
end