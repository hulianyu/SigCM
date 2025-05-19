
% obtain distance matrix using GUD distance metric

function dissim_matrix = GUD_dist( x , pm )

%% initialize the dissimilarity matrices of each attribute
dissim_matrix = cell( 1 , pm.d - pm.no_num_att ); %
for t = 1 : pm.d - pm.no_num_att %
    dissim_matrix{ t } = zeros( pm.no_values( t ) ); %
end

%% obtain conditional probability distributions 
cpd = cell( 1 , pm.d - pm.no_num_att ); %
for t = 1 : pm.d - pm.no_num_att %
    cpd{ t } = cell( pm.no_values( t ) , pm.d ); %
    for m = 1 : pm.no_values( t ) %
        locate_x_tm = x( : , t ) == m; %
        no_x_tm = sum( locate_x_tm ); %
        % for the case the related attribute is nominal and ordinal
        for r = 1 : pm.d - pm.no_num_att %
            cpd{ t }{ m , r } = zeros( 1 , pm.no_values( r ) ); %
            for g = 1 : pm.no_values( r ) %
                cpd{ t }{ m , r } ( g ) = sum( x( locate_x_tm , r ) == g ); %
            end
            cpd{ t }{ m , r } = cpd{ t }{ m , r } / no_x_tm; %
        end
        % for the case the related attribute is numerical
        for r = pm.d - pm.no_num_att + 1 : pm.d % 
            cpd{ t }{ m , r } = zeros( 1 , 5 ); %
            for s = [ 0.2, 0.4, 0.6, 0.8 ]
                cpd{ t }{ m , r }( s * 5 ) = sum( ( x( locate_x_tm , r ) >= s - 0.2 ) .* ( x( locate_x_tm , r ) < s ) );
            end
            cpd{ t }{ m , r }( 5 ) = no_x_tm - sum( cpd{ t }{ m , r } );
            cpd{ t }{ m , r } = cpd{ t }{ m , r } / no_x_tm; %
        end
    end
end

%% directly compute the cost between attributes
% initialize the interdependence matrix
inter_cost_matrix = zeros( pm.d - pm.no_num_att , pm.d );
% compute for the nominal case
for t = 1 : pm.no_nom_att
    for m = 1 : pm.no_values( t ) - 1
        for h = m + 1 : pm.no_values( t )
            cost_relate = zeros( 1 , pm.d );
            for r = 1 : pm.no_nom_att % nominal case
                diff_relate = cpd{ t }{ h , r } - cpd{ t }{ m , r };
                cost_relate( r ) = sum( abs( diff_relate ) ) / 2;
            end
            for r = pm.no_nom_att + 1 : pm.d - pm.no_num_att % ordinal case
                diff_relate = cpd{ t }{ h , r } - cpd{ t }{ m , r };
                for s = 1 : pm.no_values( r ) - 1
                    cost_relate( r ) = cost_relate( r ) + abs( diff_relate( s ) );
                    diff_relate( s +1 ) = diff_relate( s ) + diff_relate( s +1 );
                end
                cost_relate( r ) = cost_relate( r ) / ( pm.no_values( r ) - 1 );
            end
            for r = pm.d - pm.no_num_att + 1 : pm.d % numerical case
                diff_relate = cpd{ t }{ h , r } - cpd{ t }{ m , r };
                for s = 1 : 4
                    cost_relate( r ) = cost_relate( r ) + abs( diff_relate( s ) );
                    diff_relate( s +1 ) = diff_relate( s ) + diff_relate( s +1 );
                end
                cost_relate( r ) = cost_relate( r ) / 4;
            end
            inter_cost_matrix( t , : ) = inter_cost_matrix( t , : ) + cost_relate;
        end
    end
    inter_cost_matrix( t , : ) = inter_cost_matrix( t , : ) / ( ( pm.no_values( t ) - 1 ) * pm.no_values( t ) / 2 ) ;
end
% compute for the ordinal case
for t = pm.no_nom_att + 1 : pm.d - pm.no_num_att
    for m = 1 : pm.no_values( t ) - 1
        cost_relate = zeros( 1 , pm.d );
        for r = 1 : pm.no_nom_att % nominal case
            diff_relate = cpd{ t }{ m + 1 , r } - cpd{ t }{ m , r };
            cost_relate( r ) = sum( abs( diff_relate ) ) / 2;
            cost_relate( r ) = cost_relate( r ) * ( 1 / ( pm.no_values( t ) - 1 ) );
        end
        for r = pm.no_nom_att + 1 : pm.d - pm.no_num_att % ordinal case
            diff_relate = cpd{ t }{ m + 1 , r } - cpd{ t }{ m , r };
            for s = 1 : pm.no_values( r ) - 1
                cost_relate( r ) = cost_relate( r ) + abs( diff_relate( s ) );
                diff_relate( s +1 ) = diff_relate( s ) + diff_relate( s +1 );
            end
            cost_relate( r ) = cost_relate( r ) / ( pm.no_values( r ) - 1 );
            cost_relate( r ) = cost_relate( r ) * ( 1 / ( pm.no_values( t ) - 1 ) );
        end
        for r = pm.d - pm.no_num_att + 1 : pm.d % numerical case
            diff_relate = cpd{ t }{ m + 1 , r } - cpd{ t }{ m , r };
            for s = 1 : 4
                cost_relate( r ) = cost_relate( r ) + abs( diff_relate( s ) );
                diff_relate( s +1 ) = diff_relate( s ) + diff_relate( s +1 );
            end
            cost_relate( r ) = cost_relate( r ) / 4;
             cost_relate( r ) = cost_relate( r ) * ( 1 / ( pm.no_values( t ) - 1 ) );
        end
        inter_cost_matrix( t , : ) = inter_cost_matrix( t , : ) + cost_relate;
    end
end

%% directly compute the dis matrix
% for the case the target attribute is nominal
for t = 1 : pm.no_nom_att
    for m = 1 : pm.no_values( t ) - 1
        for h = m + 1 : pm.no_values( t )
            cost_relate = zeros( 1 , pm.d );
            for r = 1 : pm.no_nom_att
                diff_relate = cpd{ t }{ h , r } - cpd{ t }{ m , r };
                cost_relate( r ) = sum( abs( diff_relate ) ) / 2;
                %cost_relate( r ) = cost_relate( r ) * ( 1 / ( pm.no_values( t ) - 1 ) );
            end
            for r = pm.no_nom_att + 1 : pm.d - pm.no_num_att % ordinal case
                diff_relate = cpd{ t }{ h , r } - cpd{ t }{ m , r };
                for s = 1 : pm.no_values( r ) - 1
                    cost_relate( r ) = cost_relate( r ) + abs( diff_relate( s ) );
                    diff_relate( s +1 ) = diff_relate( s ) + diff_relate( s +1 );
                end
                cost_relate( r ) = cost_relate( r ) / ( pm.no_values( r ) - 1 );
                %cost_relate( r ) = min( ( 1 / ( pm.no_values( t ) - 1 ) ) , cost_relate( r ) );
            end
            for r = pm.d - pm.no_num_att + 1 : pm.d % numerical case
                diff_relate = cpd{ t }{ h , r } - cpd{ t }{ m , r };
                for s = 1 : 4
                    cost_relate( r ) = cost_relate( r ) + abs( diff_relate( s ) );
                    diff_relate( s +1 ) = diff_relate( s ) + diff_relate( s +1 );
                end
                cost_relate( r ) = cost_relate( r ) / 4;
            end
            dissim_matrix{ t }( m , h ) = mean( cost_relate .* inter_cost_matrix( t , : ) );
            dissim_matrix{ t }( h , m ) = dissim_matrix{ t }( m , h );
        end
    end
end
% for the case the target attribute is ordinal
for t = pm.no_nom_att + 1 : pm.d - pm.no_num_att
    dissim_vct = zeros( 1 , pm.no_values( t ) - 1 );
    for m = 1 : pm.no_values( t ) - 1
        cost_relate = zeros( 1 , pm.d );
        for r = 1 : pm.no_nom_att
            diff_relate = cpd{ t }{ m + 1 , r } - cpd{ t }{ m , r };
            cost_relate( r ) = ( sum( abs( diff_relate ) ) / 2 );
            cost_relate( r ) = cost_relate( r ) * ( 1 / ( pm.no_values( t ) - 1 ) );
        end
        for r = pm.no_nom_att + 1 : pm.d - pm.no_num_att % ordinal case
            diff_relate = cpd{ t }{ m + 1 , r } - cpd{ t }{ m , r };
            for s = 1 : pm.no_values( r ) - 1
                cost_relate( r ) = cost_relate( r ) + abs( diff_relate( s ) );
                diff_relate( s +1 ) = diff_relate( s ) + diff_relate( s +1 );
            end
            cost_relate( r ) = cost_relate( r ) / ( pm.no_values( r ) - 1 );
            cost_relate( r ) = cost_relate( r ) * ( 1 / ( pm.no_values( t ) - 1 ) );
        end
        for r = pm.d - pm.no_num_att + 1 : pm.d % numerical case
            diff_relate = cpd{ t }{ m + 1 , r } - cpd{ t }{ m , r };
            for s = 1 : 4
                cost_relate( r ) = cost_relate( r ) + abs( diff_relate( s ) );
                diff_relate( s +1 ) = diff_relate( s ) + diff_relate( s +1 );
            end
            cost_relate( r ) = cost_relate( r ) / 4;
            cost_relate( r ) = cost_relate( r ) * ( 1 / ( pm.no_values( t ) - 1 ) );
        end
        dissim_vct( m ) = mean( cost_relate .* inter_cost_matrix( t , : ) );
    end
    for m = 1 : pm.no_values( t ) - 1
        for h = m + 1 : pm.no_values( t )
            dissim_matrix{ t }( m , h ) = sum( dissim_vct( m : h - 1 ) );
            dissim_matrix{ t }( h , m ) = dissim_matrix{ t }( m , h );
        end
    end
end
for t = pm.no_nom_att + 1 : pm.d - pm.no_num_att
    dissim_matrix{ t } = dissim_matrix{ t } / max( max( ( dissim_matrix{ t } ) ) ) ;
end
end


















