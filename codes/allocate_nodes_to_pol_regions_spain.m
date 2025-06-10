%% allocate each node of the graph to a political unit

clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;

%% load country list

country_list_short;

%% choose EGM version

EGM = 8;

%% create polygon with country boundaries for each country and their union

    country_n=8; % Spain
    
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
 
    %% load political boundaries
    
    %intranational boundaries
    file = [datafolder_nuts2,'NUTS_RG_01M_2016_4326_LEVL_2.shp'];
    admin = shaperead(file,'Selector',{@(name)strcmp( name,'ES' ),'CNTR_CODE'});
     
    % keep the states within the country bounds and assign centroid to NUTS2  
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

 
%% actual bilateral trade matrix
    trade_matrix_nuts2 = csvread([datafolder_spain_trade,'trade_matrix_spain_nuts2.csv']);
    trade_matrix_nuts1 = csvread([datafolder_spain_trade,'trade_matrix_spain_nuts1.csv']);
    nuts_list = csvread([datafolder_spain_trade,'spain_nuts2.csv']);

%% represent admin as table and sort by NUTS_ID

Tadmin = struct2table(admin);
Tadmin = sortrows(Tadmin,'NUTS_ID');
Tadmin.nuts2 = nuts_list(:,1);
Tadmin.nuts1 = nuts_list(:,2);

% back to map structure
admin = table2struct(Tadmin);

% centroid of nuts1
for j=1:length(admin)    
    [admin(j).Xcent_nuts1,admin(j).Ycent_nuts1] = centroid( polyshape([admin([admin.nuts1]==admin(j).nuts1).X],[admin([admin.nuts1]==admin(j).nuts1).Y]) );    
end

%% matrixes with bilateral distances

dist_matrix.nuts2 = zeros( length(nuts_list(:,1)),length(nuts_list(:,1)) );
for i=1:length(nuts_list(:,1))
    for j=1:length(nuts_list(:,1))
        
        dist_matrix.nuts2(i,j) = deg2km( distance( admin(i).Ycent_nuts2,admin(i).Xcent_nuts2,...
                                                   admin(j).Ycent_nuts2,admin(j).Xcent_nuts2 ) );
                                               
    end
end

% list of unique nuts1 in nuts2 file
[~,IA,~] = unique([admin.nuts1]);

Xcent_nuts1 = [admin(IA).Xcent_nuts1]';
Ycent_nuts1 = [admin(IA).Ycent_nuts1]';

dist_matrix.nuts1 = zeros( length(unique(nuts_list(:,2))),length(unique(nuts_list(:,2))) );
for i=1:length(unique(nuts_list(:,2)))
    for j=1:length(unique(nuts_list(:,2)))
        dist_matrix.nuts1(i,j) = deg2km( distance( Ycent_nuts1(i),Xcent_nuts1(i),Ycent_nuts1(j),Xcent_nuts1(j) ) );
    end
end


%% load country graph

% parameters of grid to load
x_ver_hor = 0.6;   
x_diag = 0.6;      

% load spain
load( [ path_load_grids,country,'_grid_',...
                                num2str( x_diag ),'_',...
                                num2str( x_ver_hor ),'_EGM8.mat' ]  ); % this loads country_graph structure
                            
                            
%% link each node to a nuts2 and nuts1

% graph for spain
graph_export = country_graph.graph_export;

% find the nodes in each nuts2
in_nuts = zeros( graph_export.J,size( Tadmin,1 ) );
for j=1:size( Tadmin,1 )   
    in_nuts(:,j) = inpolygon( graph_export.x,graph_export.y,Tadmin.X{j},Tadmin.Y{j} );
end

% associate the nuts2 name to each node
nodes_nuts = {}; %zeros(g.J,3);
for i=1:graph_export.J
    graph_export.nodes{i}.nuts2 = nuts_list( logical( in_nuts(i,:) ),1 );
    graph_export.nodes{i}.nuts1 = nuts_list( logical( in_nuts(i,:) ),2 );
    graph_export.nodes{i}.nuts_name = convertCharsToStrings( Tadmin.NUTS_NAME{logical( in_nuts(i,:) )} );
    
    nodes_nuts{i,1} = graph_export.nodes{i}.nuts2;
    nodes_nuts{i,2} = graph_export.nodes{i}.nuts1;
    nodes_nuts{i,3} =  graph_export.nodes{i}.nuts_name;
    
    nuts1_list(i) = graph_export.nodes{i}.nuts1;
    nuts2_list(i) = graph_export.nodes{i}.nuts2;
end

graph_export.nodes_nuts = nodes_nuts;

% for each node, indicate whether neighbors are in the same nuts
for i=1:graph_export.J
        graph_export.nodes{i}.neighbors_same_nuts1 = ( nuts1_list( graph_export.nodes{i}.neighbors ) == graph_export.nodes{i}.nuts1 );
        graph_export.nodes{i}.neighbors_same_nuts2 = ( nuts2_list( graph_export.nodes{i}.neighbors ) == graph_export.nodes{i}.nuts2 );
end

% save country_graph
country_graph.graph_export = graph_export;

% save trade matrix and internal trade share
int_trade_share_nuts1 = sum( diag( trade_matrix_nuts1 ) )/...
                        sum( sum( trade_matrix_nuts1 ) );
int_trade_share_nuts2 = sum( diag( trade_matrix_nuts2 ) )/...
                        sum( sum( trade_matrix_nuts2 ) );

country_graph.trade_matrix.nuts1 = trade_matrix_nuts1;
country_graph.trade_matrix.int_trade_share_nuts1 = int_trade_share_nuts1;

country_graph.trade_matrix.nuts2 = trade_matrix_nuts2;
country_graph.trade_matrix.int_trade_share_nuts2 = int_trade_share_nuts2;

country_graph.dist_matrix = dist_matrix;

save( [ path_save_grids,country,'_grid_',...
                num2str( x_diag ),'_',...
                num2str( x_ver_hor ),'_EGM8.mat' ],'country_graph' );