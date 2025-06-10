%% create discretized maps and graphs of road network for connected countries
tic
clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;

%% code control

    % which version of EGM to use?
    EGM = 8;

    % which cities to keep?
    popthreshold = 0;         % only load places with population above threshold
    pick_largest = 0;         % =1 to use most populated place in a cell to locate the node in the cell, otherwise choose centroid
    max_N_cities = inf;       % keep the max_N_cities largest cities
    N_cat = 9;                % for plotting only. number of population quantiles is N_pop_cat+1.

    % discretization
    % these parameters control the discretization of the road network
    diagonals = 1;    % include diagonals 
    x_ver_hor = 0.6;  % tolerance for vertical-horizonal links (higher=fewer links), from 0 to 1
    x_diag = 0.6;     % tolerance for diagonal links (higher=fewer links), from 0 to 1
        
    % transform roads into graph?
    road_network_graph = 1;  % to transform roads into graphs

    % add unconnected nodes?
    N_add_edges = 0;   % number of times the original network is made denser (times unconnected nodes are connected)

    % keep just fully connected nodes?
    keep_connected = 0;

    % stock all code control variables
    code_control.datafolder_ne = datafolder_ne;
    code_control.popthreshold = popthreshold;       
    code_control.pick_largest = pick_largest;        
    code_control.max_N_cities = max_N_cities;    
    code_control.N_cat = N_cat;               
    code_control.diagonals = diagonals;
    code_control.x_ver_hor = x_ver_hor;  
    code_control.x_diag = x_diag;
    code_control.road_network_graph = road_network_graph;
    code_control.N_add_edges = N_add_edges;   
    code_control.keep_connected = keep_connected;

 
%% load country boundaries    
    file = [ datafolder_connected_countries,'connected_countries_bound_single.mat'];
    load(file);
            
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% load populated places %%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % load Population Data
   file  = [datafolder_connected_countries,'popplaces_sedac.mat'];   % load SEDAC
   temp = load(file);   
   places_sedac = temp.places_sedac;
   clear temp
     
   % total population
   totpop = sum( cell2mat( {places_sedac.population} ) );

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% bring income data %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    file  = [ datafolder_connected_countries,'places_gecon.mat'];
    load( file,'places_gecon' );
    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% bring altitude data %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    file  = [datafolder_connected_countries,'relief_etopo.mat'];
    load( file,'relief_etopo' );     
   
    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% discretize the map %%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [ code_control,places_grid,gridmap,edges,places_gecon ] = create_gridmap_connected( code_control,country_bounds,places_sedac,relief_etopo,places_gecon );
    
        % places_grid is a structure with points
            % with the same structure as places, but constructeed on averages
            % over the discrete grid. It also includes neighbors of each point.
        % gridmap is a structure with polygons
   
 %  %%%%%%%%%%%%%%%%%%%%
    %%%% load roads %%%% 
    %%%%%%%%%%%%%%%%%%%%

    % load the corresponding countries
    rdfilter = @(icc) rdfilter_connected(icc);
        
    % load the data
        load_attributes_8 = {'Shape_Leng','ICC','EXS','LTN','MED','RTT','RST','Shape_Leng'};
        file  = [ datafolder_Eurogeog_Europe_old,'RoadL.shp' ];  
        roads_temp1 = shaperead( file,'BoundingBox',country_bounds.BoundingBox,...
                                      'Attributes',load_attributes_8,...
                                      'Selector',{rdfilter,'ICC'});

    % further bound roads
    use_hull = 1;
    roads_temp2 = bound_roads( roads_temp1,country_bounds,use_hull );       
    
%
    % add weight to each segment of road based on its characteristics to compute shortest paths
    roads = assign_weights_to_roads( roads_temp2 );
    
    % save
    file = [ datafolder_connected_countries,'conn_countries_roads.mat' ];
  
    % compute some summary statistics from the road network
    
        % total KM of roads
        tot_KM_road = sum( cell2mat( {roads.totdist} ) );

        % total KM of lanes
        tot_KM_lanes = sum( cell2mat( {roads.totlanes} ) );

        % av effective lane per km or road
        av_lane_per_km = tot_KM_lanes/tot_KM_road;

    clear roads_temp1 roads_temp2;

    % transform the roads into graph
    [ road_network,unique_nodes,unique_nodes_nat ] = transform_network_to_graph( code_control,roads,[] );

        % road_network is a graph with the road network
        % unique_nodes lists the geographic coordinates of all the nodes in the network
        
    % move centroid of node in gridmap to closest node in road network within the cell
    places2 = move_cities( gridmap,places_grid,unique_nodes,unique_nodes_nat );    
    
    % discretize roads
    [ discretized_roads,graph_export,unique_edges ] = discretize_roads( code_control,gridmap,places2,road_network,edges,unique_nodes );        
    
    % add country name
    graph_export.country = {places2.country}';
    graph_export.ccode = {places2.ccode}';
    graph_export.country_n = cell2mat({places2.country_n}');
    
    % save
    country_graph.places_grid=places_grid;
    country_graph.country_bounds=country_bounds;
    country_graph.gridmap=gridmap;
    country_graph.places2=places2;
    country_graph.discretized_roads = discretized_roads;
    country_graph.graph_export = graph_export;
    country_graph.unique_edges = unique_edges;
    country_graph.edges = edges;  
    country_graph.unique_nodes = unique_nodes; 
    
        save( [ path_save_grids,'conn_countries_grid_',...
                num2str( x_diag ),'_',...
                num2str( x_ver_hor ),'_EGM8.mat' ],'country_graph' );  
        save( [ path_save_grids,'conn_countries_roads_EGM8.mat'],'roads' );   
    
toc
        