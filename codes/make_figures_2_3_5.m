%%% this code plots the country graphs created by make_country_graphs.m

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

% what to plot?
plot_baseline = 1;
plot_cfactuals = 0;

%% load country list

country_list_short;
countries1 = countries;
countries1 = [countries1;'ALL'];

%% set discretization parameters to load the data

x_ver_hor = 0.6;   % fraction of nodes in path that must go through origin-destination cell to keep kappa link
x_diag = 0.6;      % same, for diagonal paths


%% generate figures for Spain and France

for country_n = [ 10 8 ]
    
    %% load map
    country_icc = char( countries1(country_n) );  % country ICC code
    
    country = icc2name( country_icc );
    
    clc
    disp( ['Country: ',country] )
    
    load( [ path_load_grids,country,'_grid_',...
        num2str( x_diag ),'_',...
        num2str( x_ver_hor ),'_EGM8_small_cells.mat' ]  );
    
    load( [ path_load_grids,country,'_roads_EGM8_small_cells.mat'] )
    
    % unpack
    places_grid=country_graph.places_grid;
    gridmap=country_graph.gridmap;
    places2=country_graph.places2;
    discretized_roads=country_graph.discretized_roads;
    
    %% for Eurogeog
    
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
    
    %tabulate(discretized_avI)
    discretized_avI_int = round( discretized_avI );
    discretized_avI_int( ~keep_discretized ) = 0.1;
    %tabulate(discretized_avI_int)
    
    road_use_mat = cell2mat( {roads.use} )';
    road_lanes_mat = cell2mat( {roads.lanes} )';
    
    % common parameters to all graphs
    bright_line = 0.5;  % higher=brighter
    bright_nodes = 0.4; % higher=brighter
    max_bright = 9;     % higher=darker lines
    adj_width = 0.9;
    markersize = 7;
    color_growth = [ 0 1 0 ];
    color_shrink = [ 1 0 0 ];
    
    % number of categories for plotting
    N_cat = numel( unique( cell2mat( {gridmap.quantiles} ) ) )-1;
    
    %%
    if plot_baseline
        
        
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
        
    end
    
    if plot_cfactuals
        
        adj_width = 6;
        max_bright = 1.5;
        GAMMA = 0.10;
        %CFAC_TYPE = {'exp_eng','mis_eng'};
        MOBILITY = {'on','off'};
        CONGESTION = {'on'};
        
        % parameters that remain constant throughout the loop across calibrations
        alpha = 0.4;
        beta = 0.13;
        sigma = 5;
        rho = 0;
        a = 1;
        Ngoods = 10;
        nu = 1;
        
        for mobil=MOBILITY
            
            for cong = CONGESTION
                
                if strcmp(mobil,'on')
                    
                    CFAC_TYPE = {'exp_eng','mis_eng'};
                    
                elseif strcmp(mobil,'off')
                    
                    CFACT_TYPE = {'mis_eng'};
                    
                end
                
                for cfac_type = CFAC_TYPE
                    
                    for jj=1:length(GAMMA)
                        
                        gamma = GAMMA(jj);
                        
                        %% baseline parameters
                        param = init_parameters( 'a',a,'rho',rho,'alpha',alpha,'sigma',sigma,...       % preferences and technology
                            'beta',beta,'gamma',gamma,'nu',nu,'m',ones(Ngoods+1,1),...   % transpot costs
                            'K',1,...
                            'LaborMobility',char(mobil),...
                            'N',Ngoods+1,...
                            'CrossGoodCongestion',char(cong),...
                            'TolKappa',1e-4 );
                        
                        % define file to load calibration and set diary
                        filename = [ country,...
                            '_diag',num2str( x_diag ),...
                            '_hor',num2str( x_ver_hor ),...
                            '_a',num2str( param.a ),...
                            '_rho',num2str( param.rho ),...
                            '_alpha',num2str( param.alpha ),...
                            '_sigma',num2str( param.sigma ),...
                            '_beta',num2str( param.beta ),...
                            '_gamma',num2str( param.gamma ),...
                            '_nu',num2str( param.nu ),...
                            '_mobil',num2str( param.mobility ),...
                            '_cong',num2str( param.cong ),...
                            '_ngoods',num2str( Ngoods ),...
                            '_EGM8_small_cells' ];
                        
                        % load
                        load( [  path_save_cfactuals,filename,'_cfactual_',...
                            char(cfac_type),'.mat' ] )
                        
                        % recover results and g
                        g = cfactual.g;
                        results_actual = cfactual.results_actual;
                        results_cfactual = cfactual.results_cfactual;
                        
                        % GDP per capita in traded sector in calibrated model, cfactual, and data
                        y_calib = sum( results_actual.Pjn.*results_actual.Yjn,2 )./results_actual.Lj;
                        y_exp = sum( results_cfactual.Pjn.*results_cfactual.Yjn,2 )./results_cfactual.Lj;
                        y_data = g.Y./g.L;
                        
                        % total GDP in calibration and cfactual
                        PHj_calib = (1-param.alpha) / param.alpha * results_actual.PCj .* results_actual.cj ./ cfactual.param.hj;
                        GDP_calib = sum(results_actual.Pjn.*results_actual.Yjn,2) + PHj_calib.*cfactual.param.Hj;
                        
                        PHj_exp = (1-param.alpha) / param.alpha * results_cfactual.PCj .* results_cfactual.cj ./ cfactual.param.hj;
                        GDP_exp = sum(results_cfactual.Pjn.*results_cfactual.Yjn,2) + PHj_exp.*cfactual.param.Hj;
                        
                        % population in calibrated model, cfactual, and data
                        L_calib = results_actual.Lj;
                        L_exp = results_cfactual.Lj;
                        L_data = g.L;
                        
                        % total income in calibrated model, cfactual, and data
                        Y_calib = L_calib.*y_calib;
                        Y_exp = L_exp.*y_exp;
                        Y_data = L_data.*y_data;
                        
                        % total consumption in calibrated and cfactual
                        C_calib = results_actual.Cj;
                        C_exp = results_cfactual.Cj;
                        
                        % per capita
                        c_calib = C_calib./L_calib;
                        c_exp = C_exp./L_exp;
                        
                        % welfare gain
                        welfare_gain = consumption_equivalent(cfactual.param,g,c_calib,L_calib,results_cfactual.welfare)-1;
                        
                        % exports and imports in calib
                        [ X_calib,M_calib ] = recover_X_M( results_actual,g );
                        
                        % exports and imports in cfactual
                        [ X_exp,M_exp ] = recover_X_M( results_actual,g );
                        
                        % X/GDP
                        XY_calib = X_calib./Y_calib;
                        XY_exp = X_exp./Y_exp;
                        
                        % M/GDP
                        MY_calib = M_calib./Y_calib;
                        MY_exp = M_exp./Y_exp;
                        
                        % net exports over GDP
                        nx_calib = ( X_calib-M_calib )./Y_calib;
                        nx_exp = ( X_exp-M_exp )./Y_exp;
                        
                        % fundamentals
                        H_calib = cfactual.param.Hj;
                        Z_calib = max( cfactual.param.Zjn,[],2 );
                        
                        % infrastructure long
                        avI_obs = triu(g.avI);
                        avI_obs = avI_obs(:);
                        
                        avI_exp = triu(results_cfactual.Ijk);
                        avI_exp = avI_exp(:);
                        
                        % quantiles in calibrated and counterfactual model
                        quantile_L_calib = assign_quantiles( L_calib );
                        quantile_L_exp = assign_quantiles( L_exp );
                        
                        % lanes in calibrated and counterfactual model
                        I_obs = zeros( length(unique_edges),1 );
                        I_exp = zeros( length(unique_edges),1 );
                        for j=1:length(unique_edges)
                            
                            I_obs(j) = g.avI( unique_edges(j,1),unique_edges(j,2) );
                            I_exp(j) = results_cfactual.Ijk( unique_edges(j,1),unique_edges(j,2) );
                            
                        end
                        
                        % identify producers of differentiated goods
                        diff_producers = find( results_actual.Yjn(:,1)==0 );  % indexes of differentiated producers
                        
                        % compute expansion
                        gL = ( L_exp-L_calib )./( 0.5*( L_exp+L_calib ) );
                        gc = ( c_exp-c_calib )./( 0.5*( c_exp+c_calib ) );
                        gy = ( y_exp-y_calib )./( 0.5*( y_exp+y_calib ) );
                        gXY = ( XY_exp-XY_calib )./( 0.5*( XY_exp+XY_calib ) );
                        gMY = ( MY_exp-MY_calib )./( 0.5*( MY_exp+MY_calib ) );
                        gnx = ( nx_exp-nx_calib )./( 0.5*( nx_exp+nx_calib ) );
                        dI = I_exp-I_obs;
                        
                        % max I
                        maxI_obs = max( max(I_obs),max(I_exp) );
                        
                        %% plot change in network and in population
                        
                        close(figure(6))
                        h=figure(6);
                        set(h,'Position',[0 0 graph_width graph_height]);
                        
                        mapshow( country_bounds,...
                            'FaceColor','Black',...
                            'EdgeColor','White','LineWidth', 0.5 )
                        
                        % population
                        x = gL;
                        
                        % lanes
                        for j=1:length( dI )
                            
                            % COLOR SCALE OF LANES
                            colorweight = abs( dI(j) )/max(abs(dI))/max_bright;
                            lanewidth = abs( dI(j) )/max(abs(dI))*adj_width;
                            
                            if dI(j)>0
                                color =  min( bright_line+( 1-bright_line )*colorweight*color_growth,1 );
                            elseif dI(j)<0
                                color =  min( bright_line+( 1-bright_line )*colorweight*color_shrink,1 );
                            end
                            
                            mapshow( discretized_roads( j ),...
                                'Color',min( color,1 ),...
                                'LineWidth',lanewidth )
                            
                            hold on;
                            
                        end
                        
                       
                            for i=1:length( gridmap )
                                
                                norm_neg = 1/abs( mean(x(x<0)) )*1/2.5;   %gL(gL<0)*norm_neg
                                norm_pos = 1/max( x(x>0) )*2.5;         %gL(gL>0)*norm_neg
                                
                                relsize = x(i);
                                bring_fire = 0;
                                
                                if param.mobility
                                    color = min( max( -relsize*norm_neg,0 )*color_shrink+...
                                        ( max( -relsize*norm_neg,0.5 )-0.5 )*bring_fire*[ 0 1 0 ],1 )+...
                                        min( max( relsize*norm_pos,0 )*color_growth+...
                                        ( max( relsize*norm_pos,0.5 )-0.5 )*bring_fire*[ 0 1 0 ],1 );
                                else
                                    color = [ 0 0 0 ];
                                end
                                
                                mapshow( places2( i ),...
                                    'Marker','o',...
                                    'MarkerFaceColor',color,...
                                    'MarkerEdgeColor',[ 1 1 1 ],...
                                    'MarkerSize',markersize )
                                
                            end
                        
                        margin = 0.5;
                        Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
                        Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
                        set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )
                        set(gcf, 'Color', 'w');
                        
                        if save_graphs
                            export_fig([ path_final_tables,filename,'_cfactual_gL_dI_',char(cfac_type),'.eps' ], '-depsc');
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
end
