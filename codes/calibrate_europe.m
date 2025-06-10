 %%% this code runs the counterfactuals
clear all
close all
clc
tic
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

%% folder locations

folders;

%% country list

country_list_short;

%% set parameters

%calibrate constants?
calibrate_d0 = false;   % never calibrate d0 in the European case (use that from Spain)
calibrate_d1 = false;

% compute trade distance elasticity?
compute_trade_dist_elas = false;

% define the number of largest cities in each country that produce
% country-specific good
TOP_CITIES=5;

% parameters of grid to load
x_ver_hor = 0.6;
x_diag = 0.6;

% which version of EGM to use?
EGM = 8;  % must be = 8 or 10

% which level of NUTS?
NUTS = 2;  % must be = 1 or 2

% parameters that remain constant throughout the loop across calibrations
alpha = 0.4;
beta = 0.13;
sigma = 5;
rho = 0;
a = 1;
nu = 1;

% % countries
country = 'Europe';

% set parameters to be run in the loop
GAMMA = [0.10]; %[ 0.5 1 1.5 ]*beta;
MOBIL = {'off'};
CONGESTION = {'on'};

%% loop



for cong = CONGESTION
    
    for gamma = GAMMA
        
        for mobil = MOBIL
            
            % Reset the random number generator for replicability
            rng(0);
            
            % baseline parameters
            param = init_parameters( 'a',a,'rho',rho,'alpha',alpha,'sigma',sigma,...       % preferences and technology
                'beta',beta,'gamma',gamma,'nu',nu,...   % transport costs
                'K',1,...
                'LaborMobility',char(mobil),...
                'CrossGoodCongestion',char(cong),...
                'Adigator','off',...
                'TolKappa',1e-4 );
            
            % define filename to save diary, calibration and counterfactuals
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
                '_ncities',num2str( TOP_CITIES ),...
                '_EGM',num2str(EGM)];
            
            % start diary
            diary([path_save_calibrations,filename,'_calib.log']); % keep track of matlab output in log file
            clc
            
            fprintf('----------------------------------------------------------------\n');
            fprintf('CALIBRATION - Country: %s, mobil = %.1f, cong = %d, gamma = %1.2f\n\n',country,param.mobility,param.cong,param.gamma);
            
            
            %% recover graph
            
            % load map
            load( [ path_load_grids,'conn_countries','_grid_',...
                num2str( x_diag ),'_',...
                num2str( x_ver_hor ),'_EGM8.mat' ]  );
            
            % graph
            g = country_graph.graph_export;
            
            % graph_export contains information about the attributes of the discretized road network
            % graph_export.J = number of nodes
            % graph_export.nodes{n} = structure with fields specific to the node
            %    graph_export.nodes{n}.neighbors   = 1xN row with list of neighbors
            %    graph_export.nodes{n}.neighbors.x = x-coordinate
            %    graph_export.nodes{n}.neighbors.y = y-coordinate
            %    graph_export.adjacency = adjacency matrix (zero-ones)
            %    graph_export.avI = matrix with avI, equal to 1e-4 if nolink
            %    graph_export.distance = matrix with distance between connected nodes on the underlying graph
            % note: these are distances on shortest path on the
            % network if there is a link. if no link, distance is the geographic distance
            %   graph_export.L and graph_export.Y are Jx1 columns with population and income in each cell
            %   graph_export.av_alt = altitude
            %   graph_export.sd_alt = sd of altitude
            %   graph_export.rugged = ruggedness
            %   graph_export.all_distances = J*J matrix with bilateral geographic distances
            
            % we use the notation in the paper:
            %
            %       tau( Q,I ) = delta_tau * Q^beta * I^( -gamma )
            %
            
            % Fix some issues with data
            I = find( g.distance == 0 & g.adjacency==1 );
            g.distance(I) = mean(g.distance(:));
            
            %  parametrize delta_tau
            calibrate_d0 = false;
            save_calibrated_d0 = 0;
            
            d1 = 1;
            
            if param.mobility>0
                d0 = load_d0_from_file( 'calibration.dat',beta,gamma,1,param.cong,10 );
            else
                d0 = load_d0_from_file( 'calibration.dat',beta,gamma,0,param.cong,10 );
            end
            g.delta_tau = d0*( g.distance ).^d1;
            
            
            % ------------------
            % SET UP THE ECONOMY
            % ------------------
            
            % Complete and rearrange data structure
            g.ndeg = sum(reshape(tril(g.adjacency),[g.J^2,1])); % degrees of freedom
            g.Y = double(g.Y); % IPOPT only works with double precision
            g.L = double(g.L);
            g.Y = g.Y/sum(g.Y);   % normalize GDP to 1
            g.L = g.L / sum(g.L); % normalize EU population to 1
            
           
            % Correct issues with data
            I=find(g.L==0); % identify places with no population
            g.Y(I) = 1e-5; % set their production to 1e-5
            g.L(I) = 1e-5; % and their population to 1e-5
            
            % retrieve number of locations
            param.J = g.J;
            
            % Define country list and number
            [country_indices_list,~,country_short_index]=unique(g.country_n);
            ncountries=length(country_indices_list);
            
            % Define population and housing per capita
            if param.mobility==0
                param.Lj = g.L;
                param.hj=ones(g.J,1); % if not calibrated, housing is calibrated to 1 per capita
                param.Hj=param.hj.*param.Lj;
                param.omegaj=ones(g.J,1);
            elseif param.mobility==1
                param.Hj=ones(g.J,1);
            elseif param.mobility==0.5
                param.nregions=ncountries;
                g.region=country_short_index;
                param.omegar=ones(param.nregions,1);
                param.Lr=accumarray(g.region,g.L);
                param.Hj=ones(g.J,1);
            end
            
            % Define goods and productivities            
            Ngoods=ncountries+1; % one differentiated good per country+1 agricultural
            param=init_parameters('param',param,'N',Ngoods,'m',ones(Ngoods,1)); % reinit param with the right # of goods            
            
            param.Zjn = zeros( g.J,param.N );
            param.Zjn(:,1) = 1; % all cities produce the agricultural good
            
            for c=1:length(country_indices_list) % for each country
                locations=find(g.country_n==country_indices_list(c)); % get all the locations in a single country
                
                n=1;
                while n<=TOP_CITIES && ~isempty(locations) % go down the list of largest cities
                    [Y,I]=sort(g.L(locations),'descend'); % sort by population
                    param.Zjn(locations(I(1)),1+country_short_index(locations(I(1))))=1; % largest cities produce the country specific good
                    param.Zjn(locations(I(1)),1)=0; % and they don't produce the agricultural good
                    
                    % remove largest city and its neighbors from the list
                    locations=setdiff(locations,[locations(I(1)),g.nodes{locations(I(1))}.neighbors]);
                    n=n+1;
                end
            end
            
            %% calibrate       
            verbose=true;
            calibrate();
            results_actual = results; % keep the calibrated allocation for welfare comparison
            
            % Identify the delta_i_tilde's
            if param.cong==false % no cross-good congestion
                Pjkn=repmat(permute(results.Pjn,[1 3 2]),[1 g.J 1]);
                PQ=Pjkn.*results.Qjkn.^(1+param.beta);
                delta_i_tilde = g.delta_tau.*sum(PQ+permute(PQ,[2 1 3]),3)./(g.avI.^(1+param.gamma));
            else % cross-good congestion
                PCj=repmat(results.PCj,[1 g.J]);
                matm=shiftdim(repmat(param.m,[1 g.J g.J]),1);
                cost=sum(matm.*results.Qjkn.^param.nu,3).^((param.beta+1)/param.nu);
                PQ=PCj.*cost;
                delta_i_tilde = g.delta_tau.*(PQ+PQ')./(g.avI.^(1+param.gamma));
            end
            
            delta_i_tilde(~g.adjacency)=0;
            g.delta_i_tilde = delta_i_tilde;
            g.delta_i_tilde = g.delta_i_tilde / sum( g.delta_i_tilde(:).*g.avI(:) );  % normalize to meet constraint with K=1
       
            
            %% Compute also the delta_i_engineer using Collier's formula
            g.edge_ruggedness = repmat(g.rugged,[1 g.J]).^.5.*repmat(g.rugged',[g.J 1]).^.5;
            g.edge_ruggedness( g.adjacency == 0 )=0;
            
            g.delta_i_engineer = exp( -0.11*( g.distance > 50) + 0.12*log( g.edge_ruggedness ) + log( g.distance ) );
            g.delta_i_engineer( g.adjacency == 0 ) = 0;
            g.delta_i_engineer = g.delta_i_engineer / sum( g.delta_i_engineer(:).*g.avI(:) );  % normalize to meet constraint with K=1
            
            
            %% save calibrated model
            
            calibration.g = g; % note: the deltas saved under g are not normalized
            calibration.param = param;
            calibration.results_actual = results_actual;
            
            save( [ path_save_calibrations,filename,'_calib.mat' ],'calibration' );
            
            % terminate diary (saved under calibrations folder)
            diary off;
            
            
        end
    end
end

toc
