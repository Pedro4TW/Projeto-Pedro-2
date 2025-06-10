%%% this code runs the counterfactuals

clear
close all
clc
tic
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

%% folder locations

folders;

%% set parameters

% version of EGM to load
EGM = 8;

% which level of NUTS?
NUTS = 2;  % must be = 2

% how to allocate cities?
city_allocation = 'largest cities'; % 'largest cities' or 'nuts' % default is 'largest cities'

% use large cells?
cell_size = 'benchmark';  %={'benchmark','small_cells','large_cells'}

% only use large or small cells with city_allocation = 'largest cities'
if ( strcmp(cell_size,'large_cells')||strcmp(cell_size,'small_cells') ) && ~strcmp(city_allocation,'largest cities')
    error('Large or small cells only used without NUTS')
end

% on/off switch to impose an upper bound on investments
bound_investments = 1; 

% run misallocation counterfactuals?
run_misallocation_cfactual = 1;
DELTA_TYPES_MIS = {'engineer'};  % type of delta used in the misallocation counterfactual, must be 'engineer' (based on Collier)

% run expansion counterfactual?
run_expansion_cfactual = 1;
DELTA_TYPES_EXP = {'engineer','tilde'};  % type of delta used in the expansion counterfactual, can be 'engineer' or 'tilde' (using FOC's)

% run exogenous expansion?
exp_exo = 0;

% parameters of grid
x_ver_hor = 0.6;   
x_diag = 0.6;      
    
% parameters that remain constant throughout the loop across calibrations
alpha = 0.4;
beta = 0.13;
sigma = 5;
rho = 0;
a = 1;
Ngoods = 10;  % for 'largest ciies', this is the number of goods
              % for 'nuts', the number of goods is the min between Ngoods
              % and the number of NUTS regions, but the file is loaded with
              % the name Ngoods=10
nu = 1;

%% country list
switch city_allocation
    case 'largest cities'
        if strcmp(cell_size,'benchmark')
            country_list_short;   
        elseif strcmp(cell_size,'large_cells')
            country_list_short_largecells;
        elseif strcmp(cell_size,'small_cells')
            country_list_large_countries;                
        end      
    case 'nuts'
        country_list_NUTS;
end

% countries
COUNTRIES = [ 1:length(country_names(:,3)) ];  % do all countries

% set parameters to be run in the loop
GAMMA = [ ( 0.10/0.13 )*beta,( 0.13/0.10 )*beta ];
MOBIL = {'on','off'}; 
CONGESTION = {'on'};

%% loop
for gamma = GAMMA
    
    for cong = CONGESTION

        for mobil = MOBIL   

            for country_n = COUNTRIES   
                
                %% ---------------
                % COUNTERFACTUALS
                % ---------------

                % We have three different values of delta:
                % g.delta_i_engineer 
                % g.delta_i_tilde                
                % we must determine which value g.delta_i takes before running

                % NOTE:
                % all cfactual.mat have the same structure 
                % but are saved with
                % different names depending on the counterfactual
                % structure includes:
                % cfactual.g  (the delta_i used is under cfactual.g.delta_i)
                % cfactual.param
                % cfactual.results_actual
                % cfactual.results_cfactual
                % cfactual.welfare_gain
                % cfactual.results_cfactual_exo in case of expansion cfactual
                
                % baseline parameters
                param = init_parameters( 'a',a,'rho',rho,'alpha',alpha,'sigma',sigma,...       % preferences and technology
                                         'beta',beta,'gamma',gamma,'nu',nu,'m',ones(Ngoods+1,1),...   % transpot costs
                                         'K',1,...
                                         'LaborMobility',char(mobil),...
                                         'N',Ngoods+1,...
                                         'CrossGoodCongestion',char(cong),...
                                         'TolKappa',1e-4 );

                % set country  
                country_icc = char( countries( country_n ) );  % country ICC code
                country = icc2name( country_icc );
                
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
                        '_EGM',num2str(EGM)];
                    
                if strcmp(city_allocation,'nuts')
                    
                    filename = [ filename,'_nuts',num2str(NUTS) ];

                end
                
                if strcmp(cell_size,'large_cells')
                    
                    filename = [ filename,'_large_cells'];
                    
                end
                
                if strcmp(cell_size,'small_cells')
                    
                    filename = [ filename,'_small_cells'];
                    
                end
                  
                % start diary
                diary([ path_save_cfactuals,filename,'_cfactual.log']);
                clc   
                disp( ['Country: ',country] )
                
                if run_misallocation_cfactual==1||run_expansion_cfactual==1          
               

                    %% load the calibrated model
                    load( [ path_save_calibrations,filename,'_calib.mat' ] );   % this loads the structure calibration

                    g = calibration.g;
                    param = calibration.param;
                    
                    results_actual = calibration.results_actual;

                    % store
                    cfactual.g = g;
                    cfactual.param = param;
                    cfactual.results_actual = results_actual;

                    % max iterations on kappa
                    param.MAX_ITER_KAPPA = 450;

                    %% ---------------------
                    % Misallocation exercise

                    if run_misallocation_cfactual==1                    

                        % set up 
                        param.K = 1;

                        % impose upper bound if needed
                        Iu = [];
                        if bound_investments == 1
                            Iu = max(g.avI(:))*g.adjacency;
                        end                    

                        for delta_i_mis = DELTA_TYPES_MIS
                            
                                % Display stuff
                                fprintf('-------------------------\n');
                                fprintf('MISALLOCATION EXERCISE...\n\n');

                                % set delta_i
                                g.delta_i = eval( [ 'g.delta_i_',char(delta_i_mis) ] );
                                
                                % solve counterfactual
                                results_reallocation = optimal_network( param,g,g.avI,[],Iu );
                               
                                % Compute welfare gains from actual to optimal network reallocation
                                welfare_gain_reallocation = consumption_equivalent( param,g,results_actual.cj,results_actual.Lj,results_reallocation.welfare );
                                fprintf('Welfare gain from misallocation = %2.1f%%\n',(welfare_gain_reallocation-1)*100);

                                % store and save
                                cfactual.g.delta_i = g.delta_i;
                                cfactual.results_cfactual = results_reallocation;
                                cfactual.welfare_gain = welfare_gain_reallocation;                                
                                
                                if strcmp( char(delta_i_mis),'tilde' )

                                     save( [ path_save_cfactuals,filename,'_cfactual_mis_til.mat' ],'cfactual' ); 

                                 elseif strcmp( char(delta_i_mis),'engineer' )

                                     save( [ path_save_cfactuals,filename,'_cfactual_mis_eng.mat' ],'cfactual' );

                                end
                        end
                    end

                    %% ------------------------
                    % Optimal network expansion                
                    if run_expansion_cfactual==1

                        for delta_i_exp = DELTA_TYPES_EXP

                            % Display stuff
                            fprintf('-------------------------\n');
                            fprintf('EXPANSION EXERCISE: %s case...\n\n',char(delta_i_exp));
                            
                            % set network sizes
                            param.K_old = 1;                    
                            exp_ratio = 1.5;      
                            param.K = param.K_old*exp_ratio;

                            % set our choice of delta_i to run the expansion
                            g.delta_i = eval( [ 'g.delta_i_',char(delta_i_exp) ] );                        

                            % define lower bound for network expansion counterfactuals
                            Il = g.avI;                    

                            % impose upper bound if needed
                            Iu = [];
                            if bound_investments == 1
                                Iu = exp_ratio*max(g.avI(:))*g.adjacency;
                            end

                            % compute suboptimal 50% expansion everywhere & welfare gain
                            if exp_exo

                                 results_exp_exo = solve_allocation(param,g,g.avI*exp_ratio);
                                 results_exp_exo.Ijk = g.avI*exp_ratio;
                                 welfare_gain_exo = consumption_equivalent( param,g,results_actual.cj,results_actual.Lj,results_exp_exo.welfare );
                               
                            end

                            % compute optimal expansion & welfare gain
                            results_exp = optimal_network( param,g,g.avI*exp_ratio,Il,Iu );
                            welfare_gain_exp = consumption_equivalent( param,g,results_actual.cj,results_actual.Lj,results_exp.welfare );
                            fprintf('Welfare gain from %d%% expansion = %2.1f%%\n',( exp_ratio-1 )*100,( welfare_gain_exp-1 )*100);
                            
                            %% store
                            cfactual.g.delta_i = g.delta_i;  
                            cfactual.results_cfactual = results_exp;
                            cfactual.welfare_gain = welfare_gain_exp;
                            if exp_exo
                                cfactual.results_cfactual_exo = results_exp_exo;
                            end

                             % Save expansion counterfactual
                             if strcmp( char(delta_i_exp),'tilde' )

                                 save( [ path_save_cfactuals,filename,'_cfactual_exp_til.mat' ],'cfactual' ); 

                             elseif strcmp( char(delta_i_exp),'engineer' )

                                 save( [ path_save_cfactuals,filename,'_cfactual_exp_eng.mat' ],'cfactual' );

                             end

                        end

                    end    

                    toc                    

                end
                
                diary off;

           end    

        end

    end
    
end
