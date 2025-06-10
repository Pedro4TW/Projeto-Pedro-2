%% allocate each node of the country graphs to NUTS2 political subdivisions

clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

% DEFINE FOLDER LOCATIONS
folders;

%% choose type of map, countries and NUTS level
x_diag = 0.6;
x_ver_hor = 0.6;

% use large cells?
cell_size = 'benchmark';  %={'benchmark','small_cells','large_cells'}

%% load country list
switch cell_size
    case 'benchmark'
        country_list_NUTS;
    case 'large_cells'
        country_list_short_largecells_NUTS;
    case 'small_cells'
        country_list_large_countries;
end

COUNTRIES_NUTS = 1:size( country_names,1 ); 
                                        
NUTS_LEVEL = [ 1 2 3 ];

N_nuts_all = zeros( size( country_names,1 ),3 ); % table with number of nuts by country
N_nodes_all = zeros( size( country_names,1 ),1 ); % table with number of nodes by country
cellsize_all = zeros( size( country_names,1 ),1 ); % table with cellsize by country

%% create polygon with country boundaries for each country and their union
        
for country_n=COUNTRIES_NUTS
    
        % set country folder
        country_icc = char( countries(country_n) );  % country ICC code
        country = icc2name( country_icc );           % country name

        if ~ispc()   %% we are in a mac

            % folder to save the outcome
            datafolder_Eurogeog_country = [ datafolder_Eurogeog,country_icc,'/'];

        else   %% we are in a pc   

            % folder to save the outcome
            datafolder_Eurogeog_country = [ datafolder_Eurogeog,country_icc,'\' ];

        end

        clc
        disp( ['Country: ',country] )              

        % load map
        switch cell_size
            
            case 'benchmark'
                load( [ path_load_grids,country,'_grid_',...
                                            num2str( x_diag ),'_',...
                                            num2str( x_ver_hor ),'_EGM8.mat' ]  ); % this loads country_graph structure
            case 'large_cells'
                load( [ path_load_grids,country,'_grid_',...
                                            num2str( x_diag ),'_',...
                                            num2str( x_ver_hor ),'_EGM8_large_cells.mat' ]  ); % this loads country_graph structure
       
            case 'small_cells'
                load( [ path_load_grids,country,'_grid_',...
                                            num2str( x_diag ),'_',...
                                            num2str( x_ver_hor ),'_EGM8_small_cells.mat' ]  ); % this loads country_graph structure
       
        end
                                        
        % load structures. each run over nuts will add information to these
        graph_export = country_graph.graph_export;
        places_grid = country_graph.places_grid;
        gridmap = country_graph.gridmap;
        
        % nodes and cellsize by country
        N_nodes_all( country_n )= graph_export.J;
        cellsize_all( country_n ) = country_graph.code_control.cellsize;
        
        % recent nuts indexes and distance matrix for each country
        nuts = {};
        nuts_dist_matrix = {};           
                                        
        %% loop over nuts
        
        for nuts_level = NUTS_LEVEL


            %% load nuts boundaries
            eval( ['datafolder_nuts = datafolder_nuts',num2str(nuts_level)] )
            file = [datafolder_nuts,'NUTS_RG_01M_2016_4326_LEVL_',num2str(nuts_level)','.shp'];

            if strcmp(country_icc,'CHLI')
                admin = shaperead(file,'Selector',{@(name)strcmp( name,'CH' ),'CNTR_CODE'});
            else
                admin = shaperead(file,'Selector',{@(name)strcmp( name,country_icc ),'CNTR_CODE'});
            end

            %% load political boundaries

            % keep the nuts2 within the country bounds and assign centroid
            file  = [ datafolder_Eurogeog_country,'country_bounds.mat'];
            load(file);
            keep = zeros( length(admin),1 );
            for j=1:length(admin)
                  keep(j) = ~isempty( polybool( 'intersection',...
                                      cell2mat({admin(j).X}), cell2mat({admin(j).Y}),...
                                      country_bounds.X, country_bounds.Y ) );

                  [admin(j).Xcent_nuts2,admin(j).Ycent_nuts2] = centroid(polyshape(admin(j).X,admin(j).Y));
            end
            admin = admin( logical(keep) );
            
             %% Sort by NUTS_ID
             if length(admin)>1
                 Tadmin = struct2table(admin);
                 Tadmin = sortrows(Tadmin,'NUTS_ID');
                 admin = table2struct(Tadmin);
             end

            %% matrix of distances between nuts   
            N_nuts = size( admin,1 );  % number of NUTS regions
            nuts_dist_matrix{nuts_level}=zeros( N_nuts,N_nuts );
            for i=1:N_nuts
                for j=1:N_nuts

                    nuts_dist_matrix{nuts_level}(i,j)=deg2km( distance( admin(i).Ycent_nuts2,admin(i).Xcent_nuts2,...
                                                                        admin(j).Ycent_nuts2,admin(j).Xcent_nuts2 ) );
                        
                end
            end  
            
             %% link each node to a nut

             % find the nodes in each nuts
             in_nuts = zeros( graph_export.J,N_nuts );  % this will be =1 if node j is in the NUTS
             for j=1:N_nuts  
                 in_nuts(:,j) = inpolygon( graph_export.x,graph_export.y,admin(j).X,admin(j).Y );
             end
             N_nuts_all( country_n,nuts_level ) = N_nuts;
             
             % identify cases without nuts and assign closest polygon
             no_nuts = find( sum( in_nuts,2 )==0 );
             for j=1:length( no_nuts )
                 dist_temp = zeros(N_nuts,1);
                 for jj=1:N_nuts
                     temp = dsearchn( [ admin(jj).X' admin(jj).Y' ],...
                                      [ graph_export.x( no_nuts(j) ) graph_export.y( no_nuts(j) ) ] );  % closest point in the nuts
                     dist_temp(jj) = distance( [ admin(jj).X(temp) admin(jj).Y(temp) ],...
                                               [ graph_export.x( no_nuts(j) ) graph_export.y( no_nuts(j) ) ] ); % distance to the closes point in the nuts
                 end
                 in_nuts( no_nuts(j),find(dist_temp==min(dist_temp)) ) = 1; % assign the nuts
             end

            % assign corresponding NUTS to each node
            for j=1:graph_export.J

                nuts_n_temp = find(in_nuts(j,:));  % nuts is now a number from 1 to the number of nuts

                if ~isempty(nuts_n_temp)

                    nuts.codes{j,nuts_level}=nuts_n_temp;
                    nuts.names{j,nuts_level}=admin(nuts_n_temp).NUTS_NAME;
                    
                    eval( ['places_grid(j).nuts',num2str(nuts_level),'=',num2str(nuts_n_temp),';' ] )
                    eval( ['places_grid(j).nuts',num2str(nuts_level),'_name=admin(',num2str(nuts_n_temp),').NUTS_NAME;'] )

                    eval( ['gridmap(j).nuts',num2str(nuts_level),'=places_grid(j).nuts',num2str(nuts_level),';' ] )
                    eval( ['gridmap(j).nuts',num2str(nuts_level),'_name=places_grid(j).nuts',num2str(nuts_level),'_name;'] )

                end

            end          

        end
    
        % save outcomes in graph
        graph_export.nuts = nuts;
        graph_export.nuts_dist_matrix = nuts_dist_matrix;

        %% save country_graph and grids
        
        country_graph.graph_export = graph_export;
        country_graph.places_grid = places_grid;
        country_graph.gridmap = gridmap;
        
        switch cell_size
            case 'benchmark'
                 save( [ path_save_grids,country,'_grid_',...
                    num2str( x_diag ),'_',...
                    num2str( x_ver_hor ),'_EGM8.mat' ],'country_graph' );   
            case 'large_cells'
                save( [ path_save_grids,country,'_grid_',...
                    num2str( x_diag ),'_',...
                    num2str( x_ver_hor ),'_EGM8_large_cells.mat' ],'country_graph' );     
            case 'small_cells'
               save( [ path_save_grids,country,'_grid_',...
                    num2str( x_diag ),'_',...
                    num2str( x_ver_hor ),'_EGM8_small_cells.mat' ],'country_graph' );   
        end

end

%% report stats on N nuts and cell sizes for NUTS countries

        T=table(N_nuts_all,N_nodes_all,cellsize_all,country_names(:,3));
        
        switch cell_size
            case 'benchmark'
                 save( [ path_save_grids,'nuts_list.mat' ],'T' );  
            case 'large_cells'
                save( [ path_save_grids,'nuts_list_large_cells.mat' ],'T' );  
            case 'small_cells'
                save( [ path_save_grids,'nuts_list_small_cells.mat' ],'T' );  
        end

%% report nodes 
load( [ path_save_grids,'nuts_list.mat' ],'T' )