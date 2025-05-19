function result = CU(X,pi)
M = width(X);
K = max(pi); % k = 1:K
Tm = zeros(M,K);
OC = cell(M,1);
NK = sum(pi==1:K);
for m=1:M
    Xm = X(:,m);
    Q = max(Xm);
    oc = zeros(Q,K);
    for k=1:K
        set = Xm(pi==k);
        oc(:,k) = sum(set==1:Q)';  
    end
    OC{m,1} = oc;
    nq = sum(oc,2);
    Tm(m,:) = sum((oc.^2)./repmat(nq,1,K),1);
end
chi_pi = sum(Tm,1)./NK;
result = sum(chi_pi)-M;
result = result/M;
end