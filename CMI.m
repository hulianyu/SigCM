function result = CMI(X,pi)
% Cluster membership index (CMI)
N = size(X,1);
K = max(pi);
signum = zeros(K,1);
parfor k =1:K
    clusterk = find(pi==k);
    pval = SigCM_intra(X,clusterk);
    [signum(k,1),~,~] = FWER_Control(pval);
end
result = sum(signum)/N;
end