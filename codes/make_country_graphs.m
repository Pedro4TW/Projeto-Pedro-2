%% create discretized maps and graphs of road networ
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
    popthreshold = 0.6;       % only load places with population above threshold
    pick_largest = 0.6;       % =1 to use most populated place in a cell to locate the node in the cell, otherwise choose centroid
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

%% load country list

country_list_short;

%% loop through selected countries

for country_n=1:Ncountries

    clc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% load international boundaries %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % default cellsize
    country_icc = char( countries(country_n) );  % country ICC code
    country = icc2name( country_icc );           % country name

    if ~ispc()   %% we are in a mac

        % folder to save the outcome
        datafolder_Eurogeog_country =  [ datafolder_Eurogeog,country_icc,'/'];
        datafolder_Eurogeog_country_old =  [ datafolder_Eurogeog_old,country_icc,'/'];

    else   %% we are in a pc   

        % folder to save the outcome
        datafolder_Eurogeog_country =  [ datafolder_Eurogeog,country_icc,'\' ];
        datafolder_Eurogeog_country_old =  [ datafolder_Eurogeog_old,country_icc,'\' ];

    end
    
    % load country boundaries
    disp( ['Country: ',country] )
    file  = [ datafolder_Eurogeog_country,'country_bounds.mat'];
    load(file);
            
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% load populated places %%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % load Population Data
   file  = [ datafolder_Eurogeog_country,'popplaces_sedac.mat'];
   temp = load(file);   
   places_sedac = temp.places_sedac;
   clear temp
     
   % total population
   totpop = sum( cell2mat( {places_sedac.population} ) );

    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% bring income data %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    file  = [ datafolder_Eurogeog_country,'places_gecon.mat'];
    load( file,'places_gecon' );
    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% bring altitude data %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    file  = [ datafolder_Eurogeog_country,'relief_etopo.mat'];
    load( file,'relief_etopo' );     
   
    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% discretize the map %%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    default_cellsize = 0.5;   % default size is 50km by 50km
    adjust_cellsize = 1;  % if =1, scale cell size up or down to limit number of cells
    [ code_control,places_grid,gridmap,edges,places_gecon ] = create_gridmap( code_control,country_bounds,places_sedac,relief_etopo,places_gecon,default_cellsize,adjust_cellsize );
    
        % places_grid is a structure with points
            % with the same structure as places, but constructeed on averages
            % over the discrete grid. It also includes neighbors of each point.
        % gridmap is a structure with polygons
    
%%  %%%%%%%%%%%%%%%%%%%%
    %%%% load roads %%%% 
    %%%%%%%%%%%%%%%%%%%%

    % choose country to load
    icc_egm = icc2egm( country_icc );    
    rdfilter = @(icc) strcmp(icc,icc_egm);
        
    % load the data
    if EGM == 10
        load_attributes_10 = {'Shape_Leng','ICC','EXS','LTN','MED','RTT','RST','SHAPE_Leng'};
        file  = [ datafolder_Eurogeog_country,'RoadL.shp' ];  
        roads_temp1 = shaperead( file,'BoundingBox',country_bounds.BoundingBox,...
                                      'Attributes',load_attributes_10,...
                                      'Selector',{rdfilter,'ICC'});        
        [roads_temp1.('Shape_Leng')]=roads_temp1.('SHAPE_Leng');
        roads_temp1 = rmfield(roads_temp1,'SHAPE_Leng');        
    elseif EGM == 8
        load_attributes_8 = {'Shape_Leng','ICC','LTN','MED','RTT','RST','Shape_Leng'};
        file  = [ datafolder_Eurogeog_country_old,'RoadL.shp' ];  
        roads_temp1 = shaperead( file,'BoundingBox',country_bounds.BoundingBox,...
                                      'Attributes',load_attributes_8,...
                                      'Selector',{rdfilter,'ICC'});
    end

    use_hull = 0;
    roads_temp2 = bound_roads( roads_temp1,country_bounds,use_hull );    
    
    %% add weight to each segment of road based on its characteristics to compute shortest paths
    roads = assign_weights_to_roads( roads_temp2 );
    
    % compute some summary statistics from the road network
    
        % total KM of roads
        tot_KM_road = sum( cell2mat( {roads.totdist} ) );

        % total KM of lanes
        tot_KM_lanes = sum( cell2mat( {roads.totlanes} ) );

        % av effective lane per km or road
        av_lane_per_km = tot_KM_lanes/tot_KM_road;

    clear roads_temp1 roads_temp2;

    %% transform the roads into graph
    [ road_network,unique_nodes,unique_nodes_nat ] = transform_network_to_graph( code_control,roads,country );

        % road_network is the graph with the road network
        % unique_nodes lists the geographic coordinates of all the nodes in the network        
        
    %% move centroid of node in gridmap to closest node in road network within the cell
    places2 = move_cities( gridmap,places_grid,unique_nodes,unique_nodes_nat );    
    
    %% discretize roads
    [ discretized_roads,graph_export,unique_edges ] = discretize_roads( code_control,gridmap,places2,road_network,edges,unique_nodes );        
    
    %% save
    country_graph.places_grid=places_grid;
    country_graph.country_bounds=country_bounds;
    country_graph.gridmap=gridmap;
    country_graph.places2=places2;
    country_graph.discretized_roads = discretized_roads;
    country_graph.graph_export = graph_export;
    country_graph.edges = edges;  
    country_graph.unique_edges = unique_edges;      
    country_graph.unique_nodes = unique_nodes; 
    country_graph.code_control = code_control;
    
    save( [ path_save_grids,country,'_grid_',...
                num2str( x_diag ),'_',...
                num2str( x_ver_hor ),'_EGM8.mat' ],'country_graph' );   
    save( [ path_save_grids,country,'_roads_EGM8.mat'],'roads' );   
    
end
toc
        