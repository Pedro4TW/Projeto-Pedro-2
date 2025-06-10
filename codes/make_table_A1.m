%%% make table A1 in paper

clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;

% save to excel?
write_excel = 1;

%% set discretization parameters to load the data
x_ver_hor = 0.6;   % fraction of nodes in path that must go through origin-destination cell to keep kappa link
x_diag = 0.6;      % same, for diagonal paths

%% table with population, income per capita, and summary stats from road network - export to excel

    
    %% load country list
    country_list_short;
    countries1 = countries;    
    Ncountries = size(countries1,1);

    %% prepare data
    headers = { 'Country','ICC',...
                'km_of_road','Number_of_Segments','Lanes_per_km',...
                'Number_of_Cells','Total_KM_in_discretized_Network','Average_Infrastructure'};

    %% Loop
    for nn=1:Ncountries

        %% choose country
        country_icc = char( countries1(nn) );  % country ICC code    
        country = icc2name( country_icc );
        clc
        disp( ['Country: ',country] )    

        % load map
        load( [ path_save_grids,country,'_grid_',...
                num2str( x_diag ),'_',...
                num2str( x_ver_hor ),'_EGM8.mat' ] );                % this loads country_graph
     
        % load road    
        load( [ path_save_grids,country,'_roads_EGM8.mat'] );   % this loads the file roads               
        
        % unpack
        places2 = country_graph.places2;
        places_grid = country_graph.places_grid;
        discretized_roads = country_graph.discretized_roads;
        
        var =1;
        
        % country name
        COUNTRY{nn} = country;
        ICC{nn} = country_icc;
        
        % km of roads        
        KM_ACTUAL_ROADS(nn) = sum( cell2mat( {roads.totdist} ) );
        
        % line segments
        N_SEGMENTS(nn) = length(roads);
                
        % lanes per km 
        LANES_PER_KM(nn) = sum( cell2mat( {roads.totlanes} ) )/sum( cell2mat( {roads.totdist} ) );
        
        % number of cells
        N_CELLS(nn) = length( places2 );
        
        % total infrastructure
        KM_DISCRETE_ROADS(nn) = sum( cell2mat( {discretized_roads.totdist} ) );
        
        % average infrastructure
        AV_INFRASTRUCTURE(nn) = sum( cell2mat( {discretized_roads.totdist} ).*cell2mat( {discretized_roads.avI} ) )/ sum( cell2mat( {discretized_roads.totdist} ) ); 
           
        % fraction of kilometers of national roads
        length_mat =  cell2mat( {roads.totdist} );
        lanes_mat = cell2mat({roads.lanes})';
        use_mat = cell2mat({roads.use})';
        median_mat = cell2mat({roads.median})';
        paved_mat = cell2mat({roads.paved})';
        missing_use_mat = cell2mat({roads.missing_use})';
        missing_lanes_mat = cell2mat({roads.missing_use})';
        
        % average reallocation due to movement of centroid
        for j=1:length(places_grid)
            places_grid(j).X = double( places_grid(j).X );
            places_grid(j).Y = double( places_grid(j).Y );
        end
        AV_REALLOC(nn) = mean( deg2km ( distance( cell2mat({places2.Y}),cell2mat({places2.X}),...
                                     cell2mat({places_grid.Y}),cell2mat( {places_grid.X} ) ) ) );  

    end
    
%% make table and export
        
% create table
TABLE_A1 = table( COUNTRY',ICC',KM_ACTUAL_ROADS',N_SEGMENTS',LANES_PER_KM',N_CELLS',KM_DISCRETE_ROADS',AV_INFRASTRUCTURE' );

for i=1:length(headers)
    TABLE_A1.Properties.VariableNames{i} = headers{i};
end

writetable(TABLE_A1,[ path_final_tables,...
                    '/tableA1.xls'],'Sheet','raw' );