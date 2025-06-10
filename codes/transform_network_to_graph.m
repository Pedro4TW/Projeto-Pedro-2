% This code: transform shapefile with road network into a graph

function [ road_network,unique_nodes,unique_nodes_nat ] = transform_network_to_graph( code_control,roads,country )


    %% return all the coordinates defined in the road network shapefile
    all_coord  = [];  % coordinates of each node in road network
    for j=1:length( roads )

        [ j length( roads ) ]
        
        temp = [ roads(j).X( 1:length( roads(j).X )-1 )' ...
                 roads(j).Y( 1:length( roads(j).X )-1 )' ];              

        all_coord = [ all_coord;temp ];         
        
    end

    % identify the unique nodes
    [ unique_nodes,~,ic ] = unique( all_coord,'rows' );  
           % unique_nodes = coordinates of the unique nodes. nodes are numbered by its position in this vector
           % all_coord = unique_nodes( ic ), so ic indicates the unique node # associated with each elemen of all_coord
           % unique_nodes = all_coord( ia )   

    % assign an integer node index to each coordinate
    k=1;
    for j=1:length(roads)

        roads(j).node_index = ic( k:k+length( roads(j).X )-2 );   % these are the nodes indexes of all the nodes in roads(j)
        k = k+length( roads(j).X )-1;
        
    end
    
    %% return all the coordinates corresponding to national nodes defined in the road network shapefile
    
    all_coord_nat = [];  % coordinates of each node in road network    
    roads_nat = roads( logical( cell2mat( {roads.use} ) ) );   % national roads only
        
    for j=1:length( roads_nat )
        
        [ j length( roads_nat ) ]

        temp = [ roads_nat(j).X( 1:length( roads_nat(j).X )-1 )' ...
                 roads_nat(j).Y( 1:length( roads_nat(j).X )-1 )' ];              

        all_coord_nat = [ all_coord_nat;temp ];         
        
    end

    % identify the unique nodes
    [ unique_nodes_nat,~,ic_nat ] = unique( all_coord_nat,'rows' );  
           % unique_nodes = coordinates of the unique nodes. nodes are numbered by its position in this vector
           % all_coord = unique_nodes( ic ), so ic indicates the unique node # associated with each elemen of all_coord
           % unique_nodes = all_coord( ia )   

    % assign an integer node index to each coordinate
    k=1;
    for j=1:length( roads_nat )
        
        [ j length( roads_nat ) ]

        roads_nat(j).node_index = ic_nat( k:k+length( roads_nat(j).X )-2 );   % these are the nodes indexes of all the nodes in roads(j)
        k = k+length( roads_nat(j).X )-1;
        
    end

    
     
    %%
    %{
    note: 

    [ unique_nodes(roads(j).node_index,:); nan nan ] = [roads(j).X' roads(j).Y']

    

    unique_nodes( roads(j).node_index,: ) = [ roads(j).X(1:length( roads(j).X )-1 )' ...
                                              roads(j).Y(1:length( roads(j).X )-1 )']
    %}

    % define graph
    road_network_s = [];
    road_network_t = [];
    edge_lanes = [];
    edge_use = [];
    edge_distance = [];
    edge_weights = [];
    
    for j=1:length(roads)
        
        [ j length(roads) ]

        road_network_s = [ road_network_s;roads(j).node_index( 1:length( roads(j).node_index )-1 ) ];   % starting node in each edge
        
        road_network_t = [ road_network_t;roads(j).node_index( 2:length( roads(j).node_index ) ) ];     % final node in each edge
       
        edge_lanes = [ edge_lanes;ones( length( roads(j).node_index )-1,1 )*roads(j).lanes  ];  % number of actual lanes in each edge
        
        edge_use = [ edge_use;ones( length( roads(j).node_index )-1,1 )*roads(j).use  ];        % length (in KM) of each edge
        
        edge_distance = [ edge_distance;roads(j).distances' ];
        
        edge_weights = [  edge_weights;roads(j).weights' ];
    
    end


%% in countries where the same edge appears in both directions in the road data, eliminate those edges
if strcmp(country,'Finland') || strcmp(country,'Macedonia')         
    
    % repeated rows
    temp = unique( [ road_network_s road_network_t ],'rows' );    
    road_network_s = temp(:,1);
    road_network_t = temp(:,2);
    
    % repeated inverted rows
    temp = [ road_network_s road_network_t ];
    drop = [];
    for j=1:length(temp)
        for i=j+1:length(temp)
            if temp(j,1) == temp(i,2) && temp(j,2) == temp(i,1) 
                drop = [ drop;j ];           
            end
        end
    end
    
    keep = ones( length(temp),1 );
    keep(drop) = 0;
    keep = logical( keep );  
    
    road_network_s = road_network_s( keep );
    road_network_t = road_network_t( keep );    
    edge_weights = edge_weights( keep );
    edge_lanes = edge_lanes( keep );
    edge_use = edge_use( keep );
    edge_distance = edge_distance( keep );
    
end

%% graph with road network

EdgeTable = table( [road_network_s road_network_t],...   % nodes
                    edge_weights,...                     % weights
                    edge_distance,...                    % actual distance
                    edge_lanes,...              % lanes
                    edge_use,...                % use
                    'VariableNames',{'EndNodes','Weight','Distance','lanes','use'} );   % variable names
                
road_network = graph( EdgeTable );
                       
%% adjustments to the graph with road network

% connect the closest unconnected node to each node
j=0;
while j<code_control.N_add_edges
    road_network = denser_network( road_network,unique_nodes );
    j=j+1;
end

% define the network including the largest connected component
if code_control.keep_connected   % keep connected 
    keep = ( conncomp(road_network)==mode( conncomp(road_network) ) )';  % indicator for wether a node is in the largest connected set
    indexes_conn = find(keep);                                           % these are the indexes of the connected nodes in the full network
    road_network = subgraph( road_network,indexes_conn );  % note: this operation re-defines the node numbers. node numbers are different in road_network_conn than in road_network
    unique_nodes = unique_nodes( indexes_conn,: );         % coordinates of connected nodes
end