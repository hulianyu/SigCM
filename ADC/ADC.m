function pi = ADC(X, K)
% INTRODUCTION: this script conduct clustering evaluation for categorical data, i.e., data set
% contain nominal attributes only
[ pm.n , pm.d ] = size(X);
pm.no_nom_att = pm.d;
pm.no_ord_att = 0;
pm.no_num_att = 0;
pm.k = K;
for t = 1 : pm.d
    pm.no_values( t ) = length( unique( X( : , t ) ) );
end
dis_matrix = GUD_dist( X , pm );
pi = kmd_clustering( X , pm , dis_matrix );
end
