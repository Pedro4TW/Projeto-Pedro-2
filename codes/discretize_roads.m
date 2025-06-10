function [ discretized_roads,graph_export,unique_edges ] = discretize_roads( code_control,gridmap,places2,road_network,edges,unique_nodes )

%% preliminaries

% coordinates of grid centroids
Xtemp = cell2mat( {places2.X} )';
Ytemp = cell2mat( {places2.Y} )';
inroad = cell2mat( {places2.in_road_network} )';

%% vertical-horizontal neighbors

unique_edges = edges.unique_edges;
discretized_use = zeros( length( unique_edges ),1 );
discretized_dist = zeros( length( unique_edges ),1 );
discretized_lanes = zeros( length( unique_edges ),1 );
actual_dist = zeros( length( unique_edges ),1 );
keep = zeros( length( unique_edges ),1 );
paths = cell( length( unique_edges ),1 );

for j=1:length( unique_edges )
    
    [ j length( unique_edges ) ] 
    
    %% actual geographic distance in link
    actual_dist(j) = deg2km( distance( ...
                                          places2( unique_edges(j,1) ).Y,places2( unique_edges(j,1) ).X, ...
                                          places2( unique_edges(j,2) ).Y,places2( unique_edges(j,2) ).X ) );
    
    % decide whether to keep path, keep(j)=1, and return shortest path on network if so
    if inroad( unique_edges(j,1) )==1 && inroad( unique_edges(j,2) )==1  % if both are on network
        
        % origin and destination node on the road network graph
        from = places2( unique_edges(j,1) ).node;    % node in transport network corresponding to origin node in grid graph
        to = places2( unique_edges(j,2) ).node;      % node in transport network corresponding to destination node in grid graph
        
        % compute shortest path through the network
        [ path_temp,~ ] = shortestpath( road_network,from,to );
        paths{j} = path_temp;  % shortest path through actual network corresponding to unique_edges(j,:)
        
        %% plot the shortest path
        %{
             close all
             
             % show map
             mapshow( gridmap,...
            'FaceColor','Black',...
            'EdgeColor','White' ) 
        
             % plot containing cells
             mapshow( gridmap( [ unique_edges(j,1) unique_edges(j,2) ] ),...
                'FaceColor','Blue') 
            
            % plot places       
            mapshow( places2( [ unique_edges(j,1) unique_edges(j,2) ] ),...
                'Marker','o',...
                'MarkerFaceColor',[ 1 0 0 ],...
                'MarkerEdgeColor',[ 1 0 0 ],...
                'MarkerSize',3 ) 
            
            
            % show roads    
            path_temp = discretized_roads(j).opt_path;
            
            % coordinates of each link
            roads_temp.Geometry = 'Line';    
            roads_temp.X = [ unique_nodes( path_temp,1 )' NaN ];
            roads_temp.Y = [ unique_nodes( path_temp,2 )' NaN ];
            
            % plot path
            mapshow( roads_temp,...
                        'Color',...
                        discretized_roads(jj).keep*[ 0 1 0 ]+...
                        ( 1-discretized_roads(jj).keep )*[ 1 0 0 ],...
                        'LineWidth',2  );
        %}
        
        %% decide whether to keep link
        
        if ~isempty( path_temp )   % continue if there is a path

            %if code_control.keep_criterion==0

                % identify nodes belonging in source and destination cell
                nodes_in_source_cell = inpolygon( unique_nodes( path_temp,1 ),unique_nodes( path_temp,2 ),...
                                                  gridmap( unique_edges(j,1) ).X,gridmap( unique_edges(j,1) ).Y );  % nodes in path in source cell

                nodes_in_dest_cell = inpolygon( unique_nodes( path_temp,1 ),unique_nodes( path_temp,2 ),...
                                               gridmap( unique_edges(j,2) ).X,gridmap( unique_edges(j,2) ).Y ); % nodes in path in destination                            
                
                % all distances
                NN = length( path_temp ); 
                all_dist = distance( unique_nodes( path_temp( 1:NN-1 ),2 ),unique_nodes( path_temp( 1:NN-1 ),1 ),...
                                     unique_nodes( path_temp( 2:NN ),2 ),unique_nodes( path_temp( 2:NN ),1 ) );
                
                distance_in_source = sum( nodes_in_source_cell( 2:NN ).*nodes_in_source_cell( 1:NN-1 ).*all_dist )+...
                                     sum( nodes_in_source_cell( 2:NN ).*nodes_in_dest_cell( 1:NN-1 ).*all_dist );                                 
                distance_in_dest = sum( nodes_in_dest_cell( 2:NN ).*nodes_in_dest_cell( 1:NN-1 ).*all_dist )+...
                                   sum( nodes_in_dest_cell( 2:NN ).*nodes_in_source_cell( 1:NN-1 ).*all_dist );
                distance_all = sum( all_dist );                                        
                
                keep(j) = ( distance_in_source+distance_in_dest )/distance_all>=code_control.x_ver_hor; 
            
        end
        
    end

    if keep(j)==1  % link is kept
        
        % graph with the shortest path between neighboring nodes in gridmap
        shortest_path_graph = subgraph( road_network,path_temp );

        % total distance in link
        discretized_dist(j) = sum( cell2mat( table2cell( shortest_path_graph.Edges(:,'Distance') ) ) );

        % averate national use
        discretized_use(j) = ( cell2mat( table2cell( shortest_path_graph.Edges(:,'Distance') ) )'*...
                               cell2mat( table2cell( shortest_path_graph.Edges(:,'use') ) ) )/ discretized_dist(j);
                           
        % average lanes
        discretized_lanes(j) = ( cell2mat( table2cell( shortest_path_graph.Edges(:,'Distance') ) )'*...
                               cell2mat( table2cell( shortest_path_graph.Edges(:,'lanes') ) ) )/discretized_dist(j);
                           
    elseif keep(j)==0  % if no road links the two centroids
               
        % total distance in link
        discretized_dist(j) = actual_dist(j);
        
        % averate national use
        discretized_use(j) = 0;
        
        % average lanes
        discretized_lanes(j) = 0;

    end
    
end

%% diagonal neighbors

% bring in diagonal infrastructure
if code_control.diagonals == 1
    
    % diagonal distances if no direct connection through network
    unique_edges_D = edges.unique_edges_D;
    discretized_use_D = zeros( length( unique_edges_D ),1 );
    discretized_dist_D = zeros( length( unique_edges_D ),1 );
    discretized_lanes_D = zeros( length( unique_edges_D ),1 );
    actual_dist_D = zeros( length( unique_edges_D ),1 );
    keep_D = zeros( length( unique_edges_D ),1 );
    paths_D = cell( length( unique_edges_D ),1 );
    
    for j=1:length( unique_edges_D )
    
        % actual geographic distance in link
        actual_dist_D(j) = deg2km( distance( ...
                                     places2( unique_edges_D(j,1) ).Y,places2( unique_edges_D(j,1) ).X, ...
                                     places2( unique_edges_D(j,2) ).Y,places2( unique_edges_D(j,2) ).X ) );

       % decide whether to keep path, keep(j)=1, and return shortest path on network if so 
       if inroad( unique_edges_D(j,1) )==1 && inroad( unique_edges_D(j,2) )==1  % ensure that both are on network

           % origin and destination node on the road network graph
           from = places2( unique_edges_D(j,1) ).node;    % node in transport network corresponding to origin node in grid graph
           to = places2( unique_edges_D(j,2) ).node;      % node in transport network corresponding to destination node in grid graph
            
            % compute shortest path through the network
            [ path_temp,~ ] = shortestpath( road_network,from,to );         
            paths_D{j} = path_temp;
            
            %% determine whether to keep link
            
            if ~isempty( path_temp )

                %if code_control.keep_criterion==0

                    % identify nodes belonging in source and destination cell
                    nodes_in_source_cell = inpolygon( unique_nodes( path_temp,1 ),unique_nodes( path_temp,2 ),...
                                                      gridmap( unique_edges_D(j,1) ).X,gridmap( unique_edges_D(j,1) ).Y );  % nodes in path in source cell

                    nodes_in_dest_cell = inpolygon( unique_nodes( path_temp,1 ),unique_nodes( path_temp,2 ),...
                                                   gridmap( unique_edges_D(j,2) ).X,gridmap( unique_edges_D(j,2) ).Y ); % nodes in path in destination                            

                    % all distances
                    NN = length( path_temp ); 
                    all_dist = distance( unique_nodes( path_temp( 1:NN-1 ),2 ),unique_nodes( path_temp( 1:NN-1 ),1 ),...
                                         unique_nodes( path_temp( 2:NN ),2 ),unique_nodes( path_temp( 2:NN ),1 ) );

                    distance_in_source = sum( nodes_in_source_cell( 2:NN ).*nodes_in_source_cell( 1:NN-1 ).*all_dist )+...
                                         sum( nodes_in_source_cell( 2:NN ).*nodes_in_dest_cell( 1:NN-1 ).*all_dist );                                 
                    distance_in_dest = sum( nodes_in_dest_cell( 2:NN ).*nodes_in_dest_cell( 1:NN-1 ).*all_dist )+...
                                       sum( nodes_in_dest_cell( 2:NN ).*nodes_in_source_cell( 1:NN-1 ).*all_dist );
                    distance_all = sum( all_dist );
                    
                    keep_D(j) = (distance_in_source+distance_in_dest)/distance_all>=code_control.x_diag; 

            end

       end
        
       if keep_D(j)==1  % link is kept
                
                % graph with the shortest path between neighboring nodes in gridmap
                shortest_path_graph = subgraph( road_network,path_temp );

                % total distance in link
                discretized_dist_D(j) = sum( cell2mat( table2cell( shortest_path_graph.Edges(:,'Distance') ) ) );

                % averate national use
                discretized_use_D(j) = ( cell2mat( table2cell( shortest_path_graph.Edges(:,'Distance') ) )'*...
                                         cell2mat( table2cell( shortest_path_graph.Edges(:,'use') ) ) )/ discretized_dist_D(j);

                % average lanes
                discretized_lanes_D(j) = ( cell2mat( table2cell( shortest_path_graph.Edges(:,'Distance') ) )'*...
                               cell2mat( table2cell( shortest_path_graph.Edges(:,'lanes') ) ) )/discretized_dist_D(j);
                           
                                     
                                     
       elseif keep_D(j)==0  % if no road links the two centroids

                % total distance in link
                discretized_dist_D(j) = actual_dist_D(j);
                
                % averate national use
                discretized_use_D(j) = 0;
                
                % average lanes
                discretized_lanes_D(j) = 0;

       end
        
    end
    
end

%% create map structure based on discretized map

% add diagonals
if code_control.diagonals==1
    
    discretized_dist = [ discretized_dist;discretized_dist_D ];
    discretized_use = [ discretized_use;discretized_use_D ];
    discretized_lanes = [ discretized_lanes;discretized_lanes_D ];
    
    unique_edges = [ unique_edges;unique_edges_D ];
    keep = [ keep;keep_D ];
    paths = [ paths;paths_D ];
    
end

%% assign infrastructure

% average infrastructure over jk
chi_nat = 1/5.7;
avI_min = 1e-4;
discretized_avI =  discretized_lanes.*( chi_nat.^( 1-discretized_use ) );
discretized_avI( discretized_avI==0 )=avI_min;

% assign a quantile to edge
temp = discretized_avI( discretized_avI>avI_min );
quantiles_discretized_avI = quantile( temp,code_control.N_cat );

% create map structure
for j=1:length( unique_edges )
    
    % coordinates of each link
   discretized_roads(j).Geometry = 'Line';
   discretized_roads(j).X = [ Xtemp( unique_edges(j,1) ) Xtemp( unique_edges(j,2) ) NaN ];
   discretized_roads(j).Y = [ Ytemp( unique_edges(j,1) ) Ytemp( unique_edges(j,2) ) NaN ];
   
   % investment and distance in each link
   discretized_roads(j).avI = discretized_avI(j);
   discretized_roads(j).use = discretized_use(j);
   discretized_roads(j).lanes = discretized_lanes(j);
   discretized_roads(j).totdist = discretized_dist(j);
   
   % path through actual network
   discretized_roads(j).opt_path = paths{j};
   
   % is there a link?
   discretized_roads(j).keep = keep(j);
   
   % quantiles of avI of each link
   if discretized_avI(j)>avI_min
        discretized_roads(j).quantiles = sum( discretized_roads(j).avI>quantiles_discretized_avI )+1;  
   else
        discretized_roads(j).quantiles = 0.1;
   end
   
end


%% build structure to export to the solver

graph_export.J = length(places2);

x = zeros(graph_export.J,1);
y = zeros(graph_export.J,1);
for j=1:length(places2)
    
    graph_export.nodes{j}.neighbors = places2(j).neighbors; 
    graph_export.nodes{j}.x = places2(j).X;
    graph_export.nodes{j}.y = places2(j).Y;
    
    x(j) = places2(j).X;
    y(j) = places2(j).Y;
    
end

% vectors with coordinates of each location
graph_export.x=x;
graph_export.y=y;

% graph with discretized network
EdgeTable = table( unique_edges,discretized_avI,...
                   discretized_dist,discretized_use,discretized_lanes,...
                   'VariableNames',...
                   {'EndNodes','Weight',...
                   'Distance','Use','Lanes'} );   % 'Weight' is the average investment         
grid_graph = graph( EdgeTable );

% adjacency matrix - defined based on avI (variable 'Weight' defined above)
graph_export.adjacency = full( adjacency( grid_graph ) );    

% export some matrixes
nn = numnodes( grid_graph );
[s,t]=findedge( grid_graph );

    % matrix of avI
    avI_mat = full( sparse( s,t,grid_graph.Edges.Weight,nn,nn ) );  % this gives an upper triangular matrix
    avI_mat = avI_mat+avI_mat.' - diag( diag( avI_mat ) );          % this filles the lower-triangular part
    graph_export.avI = avI_mat;

    % matrix of distances
    distance_mat = full( sparse( s,t,grid_graph.Edges.Distance,nn,nn ) );       % this gives an upper triangular matrix
    distance_mat = distance_mat+distance_mat.' - diag( diag( distance_mat ) );  % this filles the lower-triangular part
    graph_export.distance = distance_mat;

%% add matrix with bilateral distances    

J = graph_export.J;
all_distances = zeros( J,J );
Xtemp = cell2mat({places2.X});
Ytemp = cell2mat({places2.Y});

for i=1:J
        all_distances( i,: ) = deg2km( distance( ones( 1,J )*Ytemp(i),ones( 1,J )*Xtemp(i),...
                                                 Ytemp,Xtemp ) );
end

% distances
graph_export.all_distances = all_distances;
    
% population 
graph_export.L = cell2mat({places2.population})';

% income
graph_export.Y = cell2mat({places2.income})';

% geography
graph_export.av_alt = cell2mat({places2.av_alt})';
graph_export.sd_alt = cell2mat({places2.sd_alt})';
graph_export.rugged = cell2mat({places2.rugged})';