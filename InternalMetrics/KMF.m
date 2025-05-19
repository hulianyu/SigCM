function result = KMF(X,pi)
result = 0;
K = length(unique(pi));
[N,M] = size(X);
for k = 1:K
    cluster_k = X(pi==k,:);
    Nk = height(cluster_k);
    mode_k = get_mode0(cluster_k);
    for i = 1:Nk
        result = result + dist_cate(cluster_k(i,:), mode_k);
    end
end
result = result/(M*N);
end

%%
function mode = get_mode0(data)
% the purpose is to obtain mode(s) of the data set, by choosing most frequent values of each column
m = length(data(1,:));
mode = zeros(1,m);
for i = 1: m
    sort_col = sort(data(:,i))';
    run_length = diff([ 0 find(sort_col(1:end-1) ~= sort_col(2:end)) length(sort_col) ]);
    uni_col = unique(sort_col);
    [~, ind] = max(run_length);
    %ind = find(run_length == val);
    mode(i) = uni_col(ind);
end
end
%%
function distance = dist_cate (a,b)
% purpose: calculate distance between two vector with distance defined by summation of
%          distance of each attribute, where distance of attribute is the indicator
%          function of non-equality.
% record of revisions:
%     date               programmer              description of change
% -----------        -----------------          ------------------------
% June 10,2003           Peng Zhang                 Original code
%
% define variables:
% a           -- vector a
% b           -- vector b
% distance    -- returned distance between a and b
%
% check whether a and b have the same length, if not set distance -1
if length(a) ~= length(b)
    distance = -1;
else
    distance = 0;
    for i = 1:length(a)
        distance = distance + (a(i) ~= b(i));
    end
end
end