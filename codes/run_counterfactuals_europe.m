%%% this code runs the counterfactuals

clear
close all
clc
tic
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;

% %% load country list
% 
% country_list_short;

%% set parameters

% version of EGM to load
EGM = 8;

% define the number of largest cities in each country that produce
% country-specific good
TOP_CITIES=5;

% on/off switch to impose an upper bound on investments
bound_investments = 1;

% run misallocation counterfactuals?
run_misallocation_cfactual = 0; %1
DELTA_TYPES_MIS = {'engineer'};  % type of delta used in the mis counterfactual, can be 'engineer', 'projected' or 'tilde'

% run expansion counterfactual?
run_expansion_cfactual = 1;
% DELTA_TYPES_EXP = {'tilde','engineer'};  % type of delta used in the exp counterfactual, can be 'engineer', 'projected' or 'tilde'
DELTA_TYPES_EXP = {'engineer'};  % type of delta used in the exp counterfactual, can be 'engineer', 'projected' or 'tilde'

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
% Ngoods = 10;
nu = 1;

% % countries
country = 'Europe';

% set parameters to be run in the loop
GAMMA = [0.10]; %[ 0.5 1 1.5 ]*beta;
MOBIL = {'partial'}; %{'off','partial','on'};
CONGESTION = {'on'};


%% loop
for cong = CONGESTION
    
    for gamma = GAMMA
        
        for mobil = MOBIL
            
            
            
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
                'beta',beta,'gamma',gamma,'nu',nu,...   % transport costs
                'K',1,...
                'LaborMobility',char(mobil),...
                'CrossGoodCongestion',char(cong),...
                'Adigator','off',...
                'TolKappa',1e-4 );
           
            
            % load calibration and set diary
            filename = [ country,...
                '_diag',num2str( x_diag ),...
                '_hor',num2str( x_ver_hor ),...
                '_a',num2str( a ),...
                '_rho',num2str( rho ),...
                '_alpha',num2str( alpha ),...
                '_sigma',num2str( sigma ),...
                '_beta',num2str( beta ),...
                '_gamma',num2str( gamma ),...
                '_nu',num2str( nu ),...
                '_mobil',num2str( param.mobility ),...
                '_cong',num2str( strcmp(cong,'on') ),...                
                '_ncities',num2str( TOP_CITIES ),...
                '_EGM',num2str(EGM) ];
            
            % start diary
            diary([ path_save_cfactuals,filename,'_cfactual.log']);
            clc
            disp( ['Country: ',country] )
            
            if run_misallocation_cfactual==1||run_expansion_cfactual==1
                
                
                %% load the calibrated model
                load( [ path_save_calibrations,filename,'_calib.mat' ] );
                
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
                        
                        % store
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
