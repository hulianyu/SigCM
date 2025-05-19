function result = Entropy(X,pi)
% expected entropy (EE) of the partition
% entropy criterion: Li T, Ma S, Ogihara M. Entropy-based criterion in categorical clustering[C]. ICML. 2004.
[N,M] = size(X);
K = max(pi);
eva_mk = zeros(M,K);
freq_pi = countcats(categorical(pi));
% eva1 = sum(times(freq_pi,log(freq_pi)));
% eva1 = eva1*M;
for m=1:M
    for k=1:K
        tmp_m = X(pi==k,m);
        freq = countcats(categorical(tmp_m));
        eva_q = times(freq,log(freq./freq_pi(k))) + times(freq_pi(k)-freq,log((freq_pi(k)-freq)./freq_pi(k)));
        tmp = sum(eva_q);
        if ~isnan(tmp)
            eva_mk(m,k) = tmp;
        end
    end
end
% eva2 = sum(eva_mk(:));
% eva = eva1 - eva2;
result = -sum(eva_mk(:))/N;
result = result/M;
end