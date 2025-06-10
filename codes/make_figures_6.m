%%% make figure 6 in the paper

clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

addpath('../Toolbox/export_fig/');

% -----------------------
% DEFINE FOLDER LOCATIONS

folders;

% save?
save_graphs = 0;

% what to plot?
plot_cfactuals = 1;

%% load country list

country_list_short;
countries1 = countries;
countries1 = [countries1;'ALL'];

%% set discretization parameters to load the data

x_ver_hor = 0.6;   % fraction of nodes in path that must go through origin-destination cell to keep kappa link
x_diag = 0.6;      % same, for diagonal paths


%% generate figure for Spain and France

country_n = 25
    
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
    discretized_ten = cell2mat( {discretized_roads.ten} )';
    discretized_lanes = cell2mat( {discretized_roads.lanes} )';
    
    max_disc_avI = max( discretized_avI );
    
    keep_discretized = cell2mat( {discretized_roads.keep} )';  % these are the links that are kept
    
    %tabulate(discretized_avI)
    discretized_avI_int = round( discretized_avI );
    discretized_avI_int( ~keep_discretized ) = 0.1;
    %tabulate(discretized_avI_int)
    
    road_use_mat = cell2mat( {roads.use} )';
    road_lanes_mat = cell2mat( {roads.lanes} )';
    road_ten_mat = cell2mat( {roads.ten} )';
    
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
    
    %%
    welfare_all = [];
    if plot_cfactuals
        
        adj_width = 1;
        GAMMA = 0.10;
        CFAC_TYPE = {'exp_eng'};   %{'exp_eng','mis_eng','exp_til'};
        MOBILITY = [ 0.5 ]; %{'on','off'};   % also: 0.5
        CONGESTION = {'on'};
        
        % parameters that remain constant throughout the loop across calibrations
        alpha = 0.4;
        beta = 0.13;
        sigma = 5;
        rho = 0;
        a = 1;
        Ngoods = 10;
        nu = 1;
        
        for cfac_type = CFAC_TYPE
            
            for cong = CONGESTION
                
                for mobil=MOBILITY
                    
                    for jj=1:length(GAMMA)
                        
                        gamma = GAMMA(jj);
                        
                        %% baseline parameters
                        
                        % define file to load calibration and set diary
                        filename = [ 'Europe',...
                            '_diag',num2str( x_diag ),...
                            '_hor',num2str( x_ver_hor ),...
                            '_a',num2str( a ),...
                            '_rho',num2str( rho ),...
                            '_alpha',num2str( alpha ),...
                            '_sigma',num2str( sigma ),...
                            '_beta',num2str( beta ),...
                            '_gamma',num2str( gamma ),...
                            '_nu',num2str( nu ),...
                            '_mobil',num2str( mobil ),...
                            '_cong',num2str( 1 ),...
                            '_ncities',num2str( 5 ),...
                            '_EGM8'];
                        
                        
                        % load
                        load( [  path_save_cfactuals,filename,'_cfactual_',...
                            char(cfac_type),'.mat' ] )
                        
                       
                        %% recover results and g
                        
                        g = cfactual.g;
                        results_actual = cfactual.results_actual;
                        results_cfactual = cfactual.results_cfactual;
                        
                        % GDP per capita in traded sector in calibrated model, cfactual, and data
                        y_calib = sum( results_actual.Pjn.*results_actual.Yjn,2 )./results_actual.Lj;
                        y_exp = sum( results_cfactual.Pjn.*results_cfactual.Yjn,2 )./results_cfactual.Lj;
                        y_data = g.Y./g.L;
                        
                        % total GDP in calibration and cfactual
                        PHj_calib = (1-alpha) / alpha * results_actual.PCj .* results_actual.cj ./ cfactual.param.hj;
                        GDP_calib = sum(results_actual.Pjn.*results_actual.Yjn,2) + PHj_calib.*cfactual.param.Hj;
                        
                        PHj_exp = (1-alpha) / alpha * results_cfactual.PCj .* results_cfactual.cj ./ cfactual.param.hj;
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
                        
                        welfare_all = [ welfare_all;welfare_gain];
                                                
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
                        
                        % check dI adds up to 0.5
                        % total_dI = sum( sum( g.delta_i.*( results_exp.I-g.avI ) ) )
                        
                        % all these cases have no construction in any case
                        % I_obs( I_obs==0.0001 )=0;
                        % I_exp( I_obs==0 )=0;
                        
                        % set maximum expansion to 10
                        %maxI_obs = max( max( I_obs ),10 );
                        %I_exp = min( I_exp,maxI_obs );  % bound counterfactual I
                        
                        % compute expansion
                        gL = ( L_exp-L_calib )./( 0.5*( L_exp+L_calib ) );
                        gc = ( c_exp-c_calib )./( 0.5*( c_exp+c_calib ) );
                        gy = ( y_exp-y_calib )./( 0.5*( y_exp+y_calib ) );
                        gXY = ( XY_exp-XY_calib )./( 0.5*( XY_exp+XY_calib ) );
                        gMY = ( MY_exp-MY_calib )./( 0.5*( MY_exp+MY_calib ) );
                        gnx = ( nx_exp-nx_calib )./( 0.5*( nx_exp+nx_calib ) );
                        dI = I_exp-I_obs;
                        dI(dI==0)=min(dI(dI>0));
                        
                        % max I
                        maxI_obs = max( max(I_obs),max(I_exp) );
                        %[ I_exp( I_obs==0.0001 ) I_obs( I_obs==0.0001 ) I_exp( I_obs==0.0001 )./I_obs( I_obs==0.0001 ) ]
                        
                        show_cells = 0;  % otherwise, show nodes
                        
                        %%
                        show_plots = 1;
                        if show_plots
                            %% plot change in network and in population
                            
                            close(figure(6))
                            h=figure(6);
                            
                            mapshow( country_bounds,...
                                'FaceColor','Black',...
                                'EdgeColor','White','LineWidth', 0.5 )
                            
                            % population
                            x = gL;
                            if show_cells
                                for i=1:length( gridmap )
                                    
                                    norm_neg = 1/abs( mean(x(x<0)) )*1/2.5;   %gL(gL<0)*norm_neg
                                    norm_pos = 1/max( x(x>0) )*2.5;         %gL(gL>0)*norm_neg
                                    
                                    relsize = x(i);
                                    bring_fire = 0;
                                    
                                    color = min( max( -relsize*norm_neg,0 )*color_shrink+...
                                        ( max( -relsize*norm_neg,0.5 )-0.5 )*bring_fire*[ 0 1 0 ],1 )+...
                                        min( max( relsize*norm_pos,0 )*color_growth+...
                                        ( max( relsize*norm_pos,0.5 )-0.5 )*bring_fire*[ 0 1 0 ],1 );
                                    
                                    mapshow( gridmap( i ),...
                                        'FaceColor',color );
                                    
                                end
                            end
                            
                            % lanes
                            for j=1:length( dI )
                                
                                % COLOR SCALE OF LANES
                                colorweight = abs( dI(j) )/max_bright;
                                
                                if dI(j)>0
                                    color =  min( bright_line+( 1-bright_line )*colorweight*color_growth,1 );
                                elseif dI(j)<0
                                    color =  min( bright_line+( 1-bright_line )*colorweight*color_shrink,1 );
                                end
                                
                                mapshow( discretized_roads( j ),...
                                    'Color',min( color,1 ),...
                                    'LineWidth',abs( dI(j) )*adj_width )
                                
                                hold on;
                                
                            end
                            
                            if ~show_cells
                                for i=1:length( gridmap )
                                    
                                    norm_neg = 1/abs( mean(x(x<0)) )*1/2.5;   %gL(gL<0)*norm_neg
                                    norm_pos = 1/max( x(x>0) )*2.5;         %gL(gL>0)*norm_neg
                                    
                                    relsize = x(i);
                                    bring_fire = 0;
                                    
                                    if mobil>0   % only show mobility in 0.5 or 1 case
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
                            end
                            
                            margin = 0.5;
                            Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
                            Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
                            set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )
                            set(gcf, 'Color', 'w');
                            
                            if save_graphs
                                
                                export_fig([  path_final_tables,filename,'_cfactual_gL_dI_',char(cfac_type),'.eps' ], '-depsc');
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end