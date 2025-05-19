% k-modes clustering algorithm
% using non-overlapped initialization
% INPUTS: data set: x
% distance/dissimilarity matrix: dis_matrix
% OUTPUTS: partition: 

function partition_new = kmd_clustering( x , pm , dis_matrix )

% initialize the k modes according to mode_initial_type
modes( 1 , : ) = x( randperm( pm.n , 1 ) , 1 : pm.d );
no_initial_mode = 1;
while no_initial_mode ~= pm.k
    no_initial_mode = no_initial_mode + 1;
    modes( no_initial_mode , : ) = x( randperm( pm.n , 1 ) , 1 : pm.d );
    if sum ( sum( modes == modes( no_initial_mode , : ) , 2 ) == pm.d ) ~= 1
        no_initial_mode = no_initial_mode - 1;
    end
end

no_total_loop = 0;
partition_old = zeros( pm.n , 1 );
obj_val = zeros( 1 , 500 );
% k-prototype clustering
loop_stop = 0;
no_loop=0;
while loop_stop == 0 && no_loop<=50
    no_total_loop = no_total_loop + 1;
    partition_new = zeros( pm.n , 1 );
    no_loop = no_loop + 1;
    dis_obj_modes = zeros( 1 , pm.k );
    for i = 1 : pm.n
        for j = 1 : pm.k
            dis_sub = zeros( 1 , pm.d );
            for t = 1 : pm.d
                dis_sub( t ) = dis_matrix{ t }( x( i , t ) , modes( j , t ) );
            end
            dis_obj_modes( j ) = sqrt(sum( dis_sub .^ 2));
        end
        [ error_sub , partition_new( i ) ] = min( dis_obj_modes );
        obj_val( no_total_loop ) = obj_val( no_total_loop ) + error_sub;      
    end
    if sum( sum( abs( partition_new - partition_old ) ) ) == 0
        loop_stop = 1;
    end
    partition_old = partition_new;
    %update modes_dct
    for c = 1 : pm.k
        x_sub = x( partition_new == c , : );
        if isempty( x_sub ) ==0
            % for nominal and ordinal attributes, cpds
            for t = 1 : pm.d
                mode_cand_dist = zeros( 1 , pm.no_values( t ) );
                modes_dct = zeros( 1 , pm.no_values( t ) );
                sta = tabulate( x_sub( : , t ) );
                sta_vct = sta( : , 3 )' / 100;
                modes_dct( 1 : length( sta_vct ) ) = sta_vct;
                for m = 1 : pm.no_values( t )
                    mode_cand_dist( m ) = ...
                        sum( dis_matrix{ t }( m , : ) .* modes_dct );
                end
                [ ~ , modes( c , t ) ] = min( mode_cand_dist );
            end
        end
    end
end
end
