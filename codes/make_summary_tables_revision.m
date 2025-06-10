%%% this code makes summary stables

clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;
warning('off')

%% tables to run
table_summary_roads = 0;
table_summary = 0;
table_all_cfactuals = 1;

% save to excel?
write_excel = 0;

%% table with summary statistic from road network
if table_summary_roads==1
    
    %% load country list
    country_list_short;
    countries1 = countries;    
    
    %% fill in data
    DATA = cell(Ncountries,26);  
    for nn=1:Ncountries

        %% choose country
        country_icc = char( countries1(nn) );  % country ICC code    
        country = icc2name( country_icc );
        clc
        disp( ['Country: ',country] )    

        % load the roads, EGM10
        temp = load( [ path_save_grids,country,'_roads.mat'] );
        roads = temp.roads;
        
        % load the roads, EGM8 
        temp = load( [ path_save_grids,country,'_roads_EGM8.mat'] );
        roads_old = temp.roads;
        
        % number of km and segments of lane by number of lanes  
        KM_lanes = zeros(1,12);
        KM_lanes_old = zeros(1,12);
        for lanes = 1:12
            select = ([roads.LTN]==lanes);
            KM_lanes(lanes) = sum([roads(select).totdist]);        
            
            select = ([roads_old.LTN]==lanes);
            KM_lanes_old(lanes) = sum([roads_old(select).totdist]);  
        end
                
        % country name
        DATA{nn,1} = country;
        DATA{nn,2} = country_icc;

        % KM of lane by lane
        for i=1:12
            DATA{nn,2+i} = KM_lanes(i);
            DATA{nn,14+i} = KM_lanes_old(i);
        end
        
    end
    
    %% write excel
    sheet_name = 'Road Network';
    xlswrite( [ path_save_summary_tables,'sumstats.xls' ],DATA,[sheet_name,'_raw'],'A1' );
    
    %% number of NUTS
    
    % benchmark
    load( [ path_save_grids,'nuts_list.mat' ],'T' );  % loads table "T"
    T_nuts = T;
    T_nuts.Properties.VariableNames{4}='Country';
    writetable( T_nuts,[ path_save_summary_tables,'sumstats.xls' ],'Sheet','number_nuts_raw','Range','A1');
    
    % large cells
    load( [ path_save_grids,'nuts_list_large_cells.mat' ],'T' );  % loads table "T"
    T_nuts = T;
    T_nuts.Properties.VariableNames{4}='Country';
    writetable( T_nuts,[ path_save_summary_tables,'sumstats.xls' ],'Sheet','number_nuts_raw','Range','H1');
    
    % small cells    
    load( [ path_save_grids,'nuts_list_small_cells.mat' ],'T' );  % loads table "T"
    T_nuts = T;
    T_nuts.Properties.VariableNames{4}='Country';
    writetable( T_nuts,[ path_save_summary_tables,'sumstats.xls' ],'Sheet','number_nuts_raw','Range','P1');
    
end       
  
%% set discretization parameters to load the data
x_ver_hor = 0.6;   % fraction of nodes in path that must go through origin-destination cell to keep kappa link
x_diag = 0.6;      % same, for diagonal paths

%% table with population, income per capita, and summary stats from road network - export to excel
if table_summary==1
    
    %% load country list
    country_list_short;
    countries1 = countries;    
    Ncountries = size(countries1,1);

    %% prepare data
    headers = { 'Country','ICC','Population','Income per Capita (thousands USD)',...
                'Number of Cells','Cell Size',...
                'km of road','Lanes per km','Average Reallocation','Number of Nodes in the Road Network','Number of Segments in the Road Network',...
                'Average Infrastructure','Total KM in discretized Network',...
                'length of national roads', 'length of non-national roads', 'av lanes of national roads', 'av lanes of non-national roads'};
    if write_excel
        xlswrite( [ path_save_summary_tables,'sumstats.xls' ],headers,'Cross-Country_raw','A1' );
    end

    countries_L = zeros( Ncountries,1 );
    countries_y = zeros( Ncountries,1 );
    countries_km_roads = zeros( Ncountries,1 );
    countries_km_lanes = zeros( Ncountries,1 );

    DATA = cell(Ncountries,10);
    shares = [];
    av_lanes = [];
    std_lanes = [];
    paved = [];
    median = [];

%% loop through countries
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
        unique_nodes = country_graph.unique_nodes;
        discretized_roads = country_graph.discretized_roads;
        
        var =1;
        
        % country name
        DATA{nn,var} = country; var=var+1;
        DATA{nn,var} = country_icc; var=var+1;

        % total population
        countries_L(nn) = sum( cell2mat( {places2.population} ) );
        DATA{nn,var} = countries_L(nn); var=var+1;

        % income per capita
        countries_y(nn) = sum( cell2mat( {places2.income} ) )/countries_L(nn);
        DATA{nn,var} = countries_y(nn)*1000000; var=var+1;

        % number of cells and cell size
        DATA{nn,var} = length( places2 ); var=var+1;
        DATA{nn,var} = country_graph.code_control.cellsize; var=var+1;
        
        % km of roads        
        countries_km_roads(nn) = sum( cell2mat( {roads.totdist} ) );
        DATA{nn,var} = countries_km_roads(nn); var=var+1;
                
        % lanes per km 
        countries_km_lanes(nn) = sum( cell2mat( {roads.totlanes} ) )/countries_km_roads(nn);
        DATA{nn,var} = countries_km_lanes(nn);    var=var+1;        

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
        DATA{nn,var} = mean( deg2km ( distance( cell2mat({places2.Y}),cell2mat({places2.X}),...
                                     cell2mat({places_grid.Y}),cell2mat( {places_grid.X} ) ) ) ); var=var+1;  
        
        % unique nodes in network
        DATA{nn,var} = length(unique_nodes); var=var+1;  
        
        % line segments
        DATA{nn,var} = length(roads); var=var+1;  
        
        % average infrastructure
        DATA{nn,var} = sum( cell2mat( {discretized_roads.totdist} ).*cell2mat( {discretized_roads.avI} ) )/ sum( cell2mat( {discretized_roads.totdist} ) ); var=var+1;  
        
        % total infrastructure
        DATA{nn,var} = sum( cell2mat( {discretized_roads.totdist} ) ); var=var+1; 
        
        % length of national roads
        DATA{nn,var} = sum( length_mat( use_mat==1 ) ); var=var+1;
        
        % length of non-national roads
        DATA{nn,var} = sum( length_mat( use_mat==0 ) ); var=var+1; 
         
        % av lanes of national roads
        DATA{nn,var} = mean( lanes_mat( use_mat==1 ) ); var=var+1; 
         
        % av lanes of non-national roads
        DATA{nn,var} = mean( lanes_mat( use_mat==0 ) ); 
                  
        % statistics of road network
        shares = [ shares; [ sum( length_mat( use_mat==0 ) ) sum( length_mat( use_mat==1 ) ) ]/sum( length_mat ) ];   
        av_lanes = [ av_lanes; [ mean( lanes_mat( use_mat==0 ) ) mean( lanes_mat( use_mat==1 ) ) ] ];
        std_lanes = [ std_lanes;[ std( lanes_mat( use_mat==0 ) ) std( lanes_mat( use_mat==1 ) ) ] ];
        paved = [ paved; [ sum( paved_mat( use_mat==0 )==1 )/length( paved_mat( use_mat==0 ) ) sum( paved_mat( use_mat==1 )==1 )/length( paved_mat( use_mat==1 ) ) ] ];
        median = [ median; [ sum( median_mat( use_mat==0 )==1 )/length( median_mat( use_mat==0 ) ) sum( median_mat( use_mat==1 )==1 )/length( median_mat( use_mat==1 ) ) ] ];
        
    end

    %% write
    if write_excel
        xlswrite( [ path_save_summary_tables,'sumstats.xls' ],DATA,'Cross-Country_raw','A2' );
    end
    
    % shares of paved/non-paved
    for i=1:2
    table_roads(:,i) = [ mean( shares(~isnan(shares(:,i)),i) );...
                   mean( av_lanes(~isnan(av_lanes(:,i)),i) );...
                   mean( std_lanes(~isnan(std_lanes(:,i)),i) );...
                   mean( median(~isnan(median(:,i)),i) );...
                   mean( paved(~isnan(paved(:,i)),i) ) ];
    end
    if write_excel
           xlswrite( [ path_save_summary_tables,'sumstats.xls' ],table_roads,'sumstats_roads_raw','A1' );
    end
    
end

%% table with all counterfactuals outcomes - export to excel

% allocate parameters that are constant over the loop

    % parameters
    alpha = 0.4;
    beta_bench = 0.13;
    sigma = 5;
    rho = 0;
    a = 1;
    nu = 1;

    % version of EGM to load
    EGM = 8;

    % which level of NUTS?
    NUTS = 2;  % must be = 1,2,3

% allocate parameters that change over the loop
    
    % beta and gamma
    GAMMA = [ ( 0.10/0.13 )*beta_bench,...
              ( 0.13/0.10 )*beta_bench ];    
    %GAMMA = [ ( 0.10/0.13 )*beta_bench ];  
          
    BETA_GAMMA = [ beta_bench GAMMA(1);...
                   beta_bench GAMMA(2) ];
  
  %  BETA_GAMMA = [ beta_bench GAMMA(1) ];


    % mobility and congestion
    MOBIL = {'on','off'};
    CONGESTION = {'on','off'};
    
    % type of cfactual
    CFAC_TYPE = {'exp_eng','mis_eng','exp_til'};
    
    % number of goods
    NGOODS = [ 10 5 15 ];

    % how to allocate cities?
    CITY_ALLOCATION = {'largest cities','nuts'}; % 'largest cities' or 'nuts'

% preallocate data
DATA = {};
cfactual_id = 1;   % unique cfactual id, common across countries
first_loop = 1;    % indicator that the first loop across countries is being run
  
%% loop baby
if table_all_cfactuals==1  

    for Ngoods = NGOODS        

        for paramset = 1:size( BETA_GAMMA,1 )

            beta = BETA_GAMMA( paramset,1 );
            gamma = BETA_GAMMA( paramset,2 );

            for cfac_type = CFAC_TYPE

                for mobil = MOBIL    

                    for cong = CONGESTION
                        
                        for city_allocation=CITY_ALLOCATION
                        
                            %% country list
                            switch char(city_allocation)
                                case 'largest cities'
                                    country_list_short;        
                                case 'nuts'
                                    country_list_NUTS;
                            end
                            countries1 = countries;    
                            Ncountries = size( countries1,1 );
                            
                            for nn = 1:Ncountries

                                 % set country  
                                 country_icc = char( countries1( nn ) );  % country ICC code
                                 country = icc2name( country_icc );
                                 clc

                                 disp( ['Country: ',country] )  

                                 param = init_parameters( 'a',a,'rho',rho,'alpha',alpha,'sigma',sigma,...       % preferences and technology
                                                                             'beta',beta,'gamma',gamma,'nu',nu,'m',ones(Ngoods+1,1),...   % transpot costs
                                                                             'K',1,...
                                                                             'LaborMobility',char(mobil),...
                                                                             'N',Ngoods+1,...
                                                                             'CrossGoodCongestion',char(cong),...
                                                                             'TolKappa',1e-4 );

                                  % define filename to load counterfactual
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
                                                '_EGM',num2str(EGM) ];
                              
                                  if strcmp( city_allocation,'nuts' )
                                      
                                      filename = [ filename,'_nuts',num2str(NUTS) ];
                                  
                                  end
                                            
                                  file_to_load = [ path_load_cfactuals,filename,'_cfactual_',...
                                                            char(cfac_type),'.mat' ];  
                                                        
                                  if isfile( file_to_load )                                                                                                            
                                      
                                      % load cfactual
                                      load( file_to_load )
                                      
                                      % recover graph and results
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

                                      % average infrastructure by location
                                      sum_avI_obs = sum( g.avI,2 );
                                      sum_avI_exp = sum( results_cfactual.Ijk,2 );

                                      % identify producers of differentiated goods
                                      diff_producers01 = ( results_actual.Yjn(:,1)==0 );  % indexes of differentiated producers = index of nodes not producing the homogeneous good

                                      % compute intra-regional trade share in counterfactual
                                      GDPj = sum(results_cfactual.Pjn.*results_cfactual.Yjn,2);
                                      intra_reg_trade = sum(GDPj);
                                      exports = sum(repmat(permute(results_cfactual.Pjn,[1 3 2]),[1 g.J 1]).*results_cfactual.Qjkn,3);
                                      results_cfactual.intra_reg_trade_share = intra_reg_trade / (intra_reg_trade + sum(exports(:)));

                                      % compute cfactual trade-distance elast
                                      results_cfactual.trade_dist_elas = trade_dist_elasticity( results_cfactual,cfactual.param,g,GDPj );
                                      
                                      % coordinates
                                      Xcoord = g.x;
                                      Ycoord = g.y;
                                      
                                      %% compute intra-regional trade share and gravity
                                      
                                      if ~max(strcmp( country_icc,{'GE','MD','ND','CY','LU','LT','LV','MK','SI'} ))   % make sure we are looking at a country that has NUTS classification
                                          
                                          % bring region and bilateral distance info
                                              if EGM == 10
                                                  load( [ path_load_grids,country,'_grid_',...
                                                      num2str( x_diag ),'_',...
                                                      num2str( x_ver_hor ),'.mat' ]  );  % load country_graph
                                              elseif EGM == 8
                                                  load( [ path_load_grids,country,'_grid_',...
                                                      num2str( x_diag ),'_',...
                                                      num2str( x_ver_hor ),'_EGM8.mat' ]  );
                                              end
                                              
                                              g.region = [ country_graph.graph_export.nuts.codes{:,NUTS} ]';  
                                              
                                              bilateral_dist = country_graph.graph_export.nuts_dist_matrix{NUTS};
                                 
                                              clear country_graph;

                                          % compute trade shares and elasticities
                                          if length(unique(g.region))==1       
                                              
                                              intra_reg_trade_share_actual = '.';
                                              intra_reg_trade_share_cfactual = '.';
                                              gravity_coeff_actual.b = '.';
                                              gravity_coeff_cfactual.b = '.';
                                          
                                          elseif length(unique(g.region))>1     
                                            
                                            [ intra_reg_trade_share_actual,gravity_coeff_actual,bilateral_trade_mat ] = compute_trade_shares_and_gravity( cfactual.param,g,results_actual,1,bilateral_dist );       % first is all, second is differentiated
                                            [ intra_reg_trade_share_cfactual,gravity_coeff_cfactual,~ ] = compute_trade_shares_and_gravity( cfactual.param,g,results_cfactual,1,bilateral_dist );                   % first is all, second is differentiated
                                                                                
                                          end
                                          
                                      else
                                          
                                          intra_reg_trade_share_actual = '.';
                                          intra_reg_trade_share_cfactual = '.';
                                          gravity_coeff_actual.b = '.';
                                          gravity_coeff_cfactual.b = '.';
                                          
                                      end
                                      
                                      
                                       %% store results

                                       if write_excel
                                           
                                           location = 1;

                                           for jjj=size(DATA,1)+1:size(DATA,1)+g.J

                                                    DATA{jjj,1} = char(country);
                                                    DATA{jjj,2} = char(country_icc);

                                                    DATA{jjj,3} = char(mobil);
                                                    DATA{jjj,4} = gamma;
                                                    DATA{jjj,5} = sigma;
                                                    DATA{jjj,6} = beta;
                                                    DATA{jjj,7} = cfactual.param.d0;
                                                    DATA{jjj,8} = cfactual.param.d1;

                                                    DATA{jjj,9} = cfactual.welfare_gain;  % note: welfare gain in % is ( welfare_gain-1 )*100)

                                                    DATA{jjj,10} = intra_reg_trade_share_actual;
                                                    DATA{jjj,11} = intra_reg_trade_share_cfactual;

                                                    DATA{jjj,12} = results_actual.trade_dist_elas;
                                                    DATA{jjj,13} = results_cfactual.trade_dist_elas;

                                                    DATA{jjj,14} = location;

                                                    DATA{jjj,15} = L_calib(location);
                                                    DATA{jjj,16} = L_exp(location);
                                                    DATA{jjj,17} = L_data(location);

                                                    DATA{jjj,18} = Y_calib(location);
                                                    DATA{jjj,19} = Y_exp(location);
                                                    DATA{jjj,20} = Y_data(location);

                                                    DATA{jjj,21} = C_calib(location);
                                                    DATA{jjj,22} = C_exp(location);

                                                    DATA{jjj,23} = XY_calib(location);
                                                    DATA{jjj,24} = XY_exp(location);
                                                    DATA{jjj,25} = MY_calib(location);
                                                    DATA{jjj,26} = MY_exp(location);

                                                    DATA{jjj,27} = H_calib(location);
                                                    DATA{jjj,28} = Z_calib(location);

                                                    DATA{jjj,29} = GDP_calib(location);
                                                    DATA{jjj,30} = GDP_exp(location);

                                                    DATA{jjj,31} = sum_avI_obs(location);
                                                    DATA{jjj,32} = sum_avI_exp(location);

                                                    DATA{jjj,33} = diff_producers01(location);

                                                    DATA{jjj,34} = g.rugged(location);

                                                    DATA{jjj,35} = char(cfac_type);

                                                    DATA{jjj,36} = char(cong);

                                                    DATA{jjj,37} = Ngoods;

                                                    DATA{jjj,38} = cfactual_id;

                                                    DATA{jjj,39} = paramset;

                                                    if isfield(results_actual,'taxes')
                                                            
                                                        DATA{jjj,40} = double(results_actual.taxes.mean_ad_valorem);

                                                        DATA{jjj,41} = double(results_actual.taxes.total_t_over_GDP);

                                                        DATA{jjj,42} = double(results_actual.taxes.total_t_over_tradeable);
                                                        
                                                    else
                                                        
                                                        DATA{jjj,40} = '.';
                                                        DATA{jjj,41} = '.';
                                                        DATA{jjj,42} = '.';
                                                        
                                                    end
                                                    
                                                    if isfield(results_actual,'flaw_coeffs')

                                                        DATA{jjj,43} = double(results_actual.flaw_coeffs.local_elasticity_Q_to_I);

                                                        DATA{jjj,44} = double(results_actual.flaw_coeffs.agg_elasticity_Q_to_I);
                                                        
                                                    else
                                                        
                                                        DATA{jjj,43} = '.';
                                                        DATA{jjj,44} = '.';
                                                        
                                                    end
                                                    
                                                    DATA{jjj,45} = char(city_allocation);
                                                    
                                                    DATA{jjj,46} = gravity_coeff_actual.b;
                                                    
                                                    DATA{jjj,47} = gravity_coeff_cfactual.b;

                                                    DATA{jjj,48} = double( Xcoord(location) );
                                                    
                                                    DATA{jjj,49} = double( Ycoord(location) );
                                                    
                                                    location = location+1;
                                           end
                                                
                                       end

                                  end
                                  
                            end

                            first_loop = 0;
                            cfactual_id = cfactual_id + 1;
                        
                        end

                    end

                end    

            end

        end

    end
    
    if write_excel

              headers =  { 'Country','ICC',...
                         'Mobility','gamma','sigma','beta','d0','d1',...
                         'welfare_gain',...
                         'intra_reg_calib','intra_reg_exp',...
                         'dist_elast_calib','dist_elast_exp',...
                         'location',...
                         'L_calib','L_exp','L_data',...
                         'Y_calib','Y_exp','Y_data',...
                         'C_calib','C_exp',...
                         'XY_calib','XY_exp','MY_calib','MY_exp',...
                         'H_calib','Z_calib',...
                         'GDP_calib','GDP_exp',...
                         'I_obs','I_exp',...
                         'diff_producer',...
                         'rugged',...
                         'cfac_type',...
                         'cong',...
                         'Ngoods',...
                         'cfactual_id',...
                         'paramset',...
                         'mean_t','total_t_over_GDP','total_t_over_T',...
                         'flaw_loc_elast','flaw_agg_elast',...
                         'city_allocation',...
                         'gravity_coeff_calib','gravity_coeff_exp',...
                         'Xcoord','Ycoord'};

                    xlswrite( [ path_save_summary_tables,'all_cfactuals.xls' ],headers,'all_results','A1' );
                    xlswrite( [ path_save_summary_tables,'all_cfactuals.xls' ],DATA,'all_results','A2' );          

       end 

end


        
        