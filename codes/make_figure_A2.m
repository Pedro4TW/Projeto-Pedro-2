% generate figure A2 in the paper

clear
close all
clc

set(0,'DefaultFigureWindowStyle','normal','DefaultFigureVisible','on') 

addpath('../Toolbox/export_fig/');

graph_width=800;
graph_height=600;

% -----------------------
% DEFINE FOLDER LOCATIONS

folders;

% save?
save_graphs=0;

%% load country list

country_list_short;
countries1 = countries;
countries1 = [countries1;'ALL'];

%% set discretization parameters to load the data

x_ver_hor = 0.6;   % fraction of nodes in path that must go through origin-destination cell to keep kappa link
x_diag = 0.6;      % same, for diagonal paths


%% generate figures for Europe

country_n = 25;
    
    %% load map
    country_icc = char( countries1(country_n) );  % country ICC code
    
    country = icc2name( country_icc );
    
    clc
    disp( ['Country: ',country] )
    
    load( [ path_load_grids,country,'_grid_',...
        num2str( x_diag ),'_',...
        num2str( x_ver_hor ),'_tentec.mat' ]  );
    
    load( [ path_load_grids,country,'_roads_tentec.mat'] )
    
    % unpack
    places_grid=country_graph.places_grid;
    gridmap=country_graph.gridmap;
    places2=country_graph.places2;
    discretized_roads=country_graph.discretized_roads;
    graph_export=country_graph.graph_export;
    unique_edges=country_graph.unique_edges;
    edges = country_graph.edges;
    unique_nodes=country_graph.unique_nodes;
    
    % country bounds
    country_bounds=country_graph.country_bounds;
    
    %% compute some variables used in figures
    
    discretized_avI = cell2mat( {discretized_roads.avI} )';
    discretized_use = cell2mat( {discretized_roads.use} )';
    discretized_lanes = cell2mat( {discretized_roads.lanes} )';
    
    max_disc_avI = max( discretized_avI );
    
    keep_discretized = cell2mat( {discretized_roads.keep} )';  % these are the links that are kept
    
    discretized_avI_int = round( discretized_avI );
    discretized_avI_int( ~keep_discretized ) = 0.1;
    
    road_use_mat = cell2mat( {roads.use} )';
    road_lanes_mat = cell2mat( {roads.lanes} )';
    
    % common parameters to all graphs
    bright_line = 0.5;  % higher=brighter
    bright_nodes = 0.4; % higher=brighter
    max_bright = 9;    % higher=darker lines
    adj_width = 0.9;
    markersize = 7;
    color_growth = [ 0 1 0 ];   % red = [ 1 0 0 ]  green = [ 0 1 0 ] blue = [ 0 0 1 ]
    color_shrink = [ 1 0 0 ];
    
    % number of categories for plotting
    N_cat = numel( unique( cell2mat( {gridmap.quantiles} ) ) )-1;
    
        
        %% population
        
        close(figure(1))
        fig=figure(1);
        set(fig,'Position',[0 0 graph_width graph_height]);
        
        % population
        for i=N_cat+1:-1:1
            
            relsize = bright_nodes+(1-bright_nodes)*i/( N_cat+1 );
            mapshow( gridmap( cell2mat( {gridmap.quantiles} )==i ),...
                'FaceColor',relsize*[ 0 1 1 ] )
        end
        mapshow( gridmap( cell2mat( {gridmap.quantiles} )==0 ),...
            'FaceColor',[ 1 1 1 ] )
        
        % centroids
        for i=N_cat+1:-1:1
            
            relsize = bright_nodes+(1-bright_nodes)*i/( N_cat+1 );
            
            mapshow( places2( cell2mat( {gridmap.quantiles} )==i ),...
                'Marker','o',...
                'MarkerFaceColor',relsize*[ 0 1 1 ],...
                'MarkerEdgeColor',[ 0 0 0 ],...
                'MarkerSize',6 )
        end
        
        margin = 0.5;
        Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
        Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
        set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )
        
        if save_graphs
            export_fig([ path_final_tables,country,'_pop.eps' ], '-depsc');
        end
        
        %% baseline grid
        close(figure(2))
        h=figure(2);
        set(h,'Position',[0 0 graph_width graph_height]);
        
        mapshow( gridmap,...
            'FaceColor','Black',...
            'EdgeColor','White','LineWidth', 0.5 )
        
        mapshow( discretized_roads,...
            'Color',[ 1 0 0 ],...
            'LineWidth', 1  )
        
        for i=N_cat+1:-1:1
            
            relsize = bright_nodes+(1-bright_nodes)*i/( N_cat+1 );
            
            mapshow( places2( cell2mat( {gridmap.quantiles} )==i ),...
                'Marker','o',...
                'MarkerFaceColor',relsize*[ 0 1 1 ],...
                'MarkerEdgeColor',[ 0 0 0 ],...
                'MarkerSize',markersize )
        end
        
        
        margin = 0.5;
        Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
        Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
        set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )
        if save_graphs
            export_fig([ path_final_tables,country,'_baseline_grid.eps' ], '-depsc');
        end
        
        %% actual centroids and actual (non-discretized) infrastructure
        
        close(figure(3))
        h=figure(3);
        set(h,'Position',[0 0 graph_width graph_height]);
        
        % show bounds
        mapshow( country_bounds,...
            'FaceColor','Black',...
            'EdgeColor','White','LineWidth', 0.5 )
        
        % show roads
        for j=1:max( road_lanes_mat )
            
            keep = logical( ( road_lanes_mat==j ).*( road_use_mat==0 ) );
            mapshow( roads( keep ),...
                'Color',min( bright_line+( 1-bright_line )*j/max_bright*[ 1 0 0 ],1 ),...
                'LineWidth', j*adj_width  );
            
            keep = logical( ( road_lanes_mat==j ).*( road_use_mat==1 ) );
            mapshow( roads( keep ),...
                'Color',min( bright_line+( 1-bright_line )*j/max_bright*[ 0 1 0 ],1 ),...
                'LineWidth', j*adj_width   );
            
        end
        
        % centroids
        for i=N_cat+1:-1:1
            
            relsize = bright_nodes+(1-bright_nodes)*i/( N_cat+1 );
            
            mapshow( places2( cell2mat( {gridmap.quantiles} )==i ),...
                'Marker','o',...
                'MarkerFaceColor',relsize*[ 0 1 1 ],...
                'MarkerEdgeColor',[ 0 0 0 ],...
                'MarkerSize',markersize )
        end
        
        margin = 0.5;
        Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
        Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
        set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )
        
        if save_graphs
            export_fig([ path_final_tables,country,'_actual_network.eps' ], '-depsc');
        end
        
        %% actual centroids and actual (discretized) infrastructure
        
        close(figure(4))
        h=figure(4);
        set(h,'Position',[0 0 graph_width graph_height]);
        
        % bounds
        mapshow( country_bounds,...
            'FaceColor','Black',...
            'EdgeColor','White','LineWidth', 0.5 )
        
        
        % show lanes and use
        for j=1:length( discretized_roads )
            
            if discretized_roads(j).keep>0
                
                colorweight = discretized_lanes(j)/max_bright;
                
                color =  bright_line + ( 1-bright_line ) * colorweight * ...
                    ( discretized_use(j)*[ 0 1 0 ] +...
                    + ( 1-discretized_use(j) )*[ 1 0 0 ] );
                
                mapshow( discretized_roads( j ),...
                    'Color',min( color,1 ),...
                    'LineWidth',discretized_lanes(j)*adj_width )
                
                hold on;
                
            end
            
        end
        
        for i=N_cat+1:-1:1
            
            relsize = bright_nodes+(1-bright_nodes)*i/( N_cat+1 );
            
            mapshow( places2( cell2mat( {gridmap.quantiles} )==i ),...
                'Marker','o',...
                'MarkerFaceColor',relsize*[ 0 1 1 ],...
                'MarkerEdgeColor',[ 0 0 0 ],...
                'MarkerSize',markersize )
        end
        
        margin = 0.5;
        Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
        Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
        set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )
        if save_graphs
            export_fig([ path_final_tables,country,'_discrete_network.eps' ], '-depsc');
        end
        
   