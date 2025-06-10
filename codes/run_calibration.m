%%% this code runs the counterfactuals

clear
close all
clc
tic
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

%% folder locations

folders;

%% set parameters

% which version of EGM to use?
EGM = 8;  % must be = 8 or 10

% How to allocate cities?
% - largest_cities: only the largest city in a country get assigned a
% differentiated good to produce.
% - nuts: there is one NUTS-specific product
city_allocation = 'largest cities'; % 'largest cities' or 'nuts'
if ~validatestring(city_allocation,{'largest cities','nuts'})
    fprintf('%s.m: unknown city allocation method ''%s.''\n',mfilename,city_allocation);
end

% which level of NUTS? (if nuts is selected in city_allocation)
NUTS = 2;  % must be = 2

% use large cells?
cell_size = 'benchmark';  %={'benchmark','small_cells','large_cells'}

% only use large or small cells with city_allocation = 'largest cities'
if ( strcmp(cell_size,'large_cells')||strcmp(cell_size,'small_cells') ) && ~strcmp(city_allocation,'largest cities')
    error('Large or small cells only used without NUTS')
end

% determine file to load d0
d0_file='calibration';
switch city_allocation
    case 'largest cities'
        if strcmp(cell_size,'large_cells')
            d0_file=[d0_file,'_large_cells'];
        elseif strcmp(cell_size,'small_cells')
            d0_file=[d0_file,'_small_cells'];
        end
    case 'nuts'
        d0_file=[d0_file,sprintf('_NUTS%d.dat',NUTS)];
end
d0_file=[d0_file,'.dat'];

%calibrate constants?
calibrate_d0 = true;      % if true, only calibrated if the country is Spain
calibrate_d1 = false;

% parameters of grid to load
x_ver_hor = 0.6;
x_diag = 0.6;

% parameters that remain constant throughout the loop across calibrations
alpha = 0.4;
beta = 0.13;
sigma = 5;
rho = 0;
a = 1;
nu = 1;

% number of differentiated goods
Ngoods = 10;  % for 'largest cities', this is the number of goods
% for 'nuts', the number of goods is the min between Ngoods and the number of NUTS regions

Ngoods_fixed = Ngoods; % this will not change throughout the code. used to define the filename and d0

% use adigator?
adigator = 0;

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

% countries for loop
n_Spain = find( strcmp([country_names(:,3)],'Spain') );
COUNTRIES = [ n_Spain 1:(n_Spain-1) (n_Spain+1):length(country_names(:,3)) ]; % start with Spain, then loop through the rest

% set parameters to be run in the loop
GAMMA = [ ( 0.10/0.13 )*beta,( 0.13/0.10 )*beta ];
MOBIL = {'on','off'};
CONGESTION = {'on'};

%% loop

for country_n=COUNTRIES
    
    for cong = CONGESTION
        
        for gamma = GAMMA
            
            for mobil = MOBIL
                
                % baseline parameters
                % this is used to create filename
                % number of goods may be redefined below to determine Ngoods in 'nuts' case
                
                param = init_parameters( 'a',a,'rho',rho,'alpha',alpha,'sigma',sigma,...
                    'beta',beta,'gamma',gamma,'nu',nu,'m',ones(Ngoods+1,1),...
                    'K',1,...
                    'LaborMobility',char(mobil),...
                    'N',Ngoods+1,...
                    'CrossGoodCongestion',char(cong),...
                    'TolKappa',1e-4, 'Adigator','off' );
                
                % set country
                country_icc = char( countries( country_n ) );  % country ICC code
                country = icc2name( country_icc );
                
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
                    '_ngoods',num2str( Ngoods_fixed ),...
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
                diary( [path_save_calibrations,filename,'_calib.log'] ); % keep track of matlab output in log file
                clc
                
                fprintf('----------------------------------------------------------------\n');
                fprintf('CALIBRATION - Country: %s, mobil = %d, cong = %d, gamma = %1.2f\n\n',country,param.mobility,param.cong,param.gamma);
                
                % calibrate d0 only if country is Spain and we run the calibration
                if strcmp(country,'Spain') && (calibrate_d0==true)
                    calibrate_d0 = true;
                    save_calibrated_d0 = 1;
                else
                    calibrate_d0 = false;
                    save_calibrated_d0 = 0;
                end
                
                %% recover graph
                
                graph_filename = [ path_load_grids,country,'_grid_',...
                    num2str( x_diag ),'_',...
                    num2str( x_ver_hor ) ];
                if EGM == 8
                    graph_filename = [ graph_filename,'_EGM8' ];
                end
                
                if strcmp(cell_size,'large_cells')
                    graph_filename = [ graph_filename,'_large_cells' ];
                end
                
                if strcmp(cell_size,'small_cells')
                    graph_filename = [ graph_filename,'_small_cells' ];
                end
                
                load( [ graph_filename,'.mat' ] );
                
                % graph
                g = country_graph.graph_export;
                unique_edges = country_graph.unique_edges;
                
                % prepare NUTS data
                if isfield(g,'nuts')
                    g.region = [g.nuts.codes{:,NUTS}]';
                    NUTS_list = unique(g.region);
                    nb_NUTS = length(NUTS_list); % # of NUTS regions
                end
                
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
                %d0 = 1; % default value of d0
                d1 = 1;
                if calibrate_d0==false
                    d0 = load_d0_from_file( d0_file,beta,gamma,param.mobility,param.cong,Ngoods_fixed );
                    g.delta_tau = d0*( g.distance ).^d1;
                end
                
                % ------------------
                % SET UP THE ECONOMY
                % ------------------
                
                % Complete and rearrange data structure
                g.ndeg = sum(reshape(tril(g.adjacency),[g.J^2,1])); % degrees of freedom
                g.Y = double(g.Y); % IPOPT only works with double precision
                g.L = double(g.L);
                g.Y = g.Y/sum(g.Y);   % normalize GDP to 1
                g.L = g.L/sum(g.L); % normalize population to 1
                
                % Correct issues with data
                I=find(g.L==0); % identify places with no population
                g.Y(I) = 1e-5;  % set their production to 1e-5
                g.L(I) = 1e-5;  % and their population to 1e-5
                
                % retrieve number of locations
                param.J = g.J;
                
                % Define population and housing per capita
                if param.mobility==false
                    param.Lj = g.L;
                    param.hj=ones(g.J,1); % if not calibrated, housing is calibrated to 1 per capita
                    param.Hj=param.hj.*param.Lj;
                    param.omegaj=ones(g.J,1);
                else
                    param.Hj=ones(g.J,1);
                end
                
                % Define goods and productivities
                switch city_allocation
                    
                    case 'largest cities'
                        
                        param.Zjn = zeros( g.J,Ngoods+1 );
                        param.Zjn(:,1) = 1; % all cities produce the agricultural good
                        
                        locations=1:g.J;
                        n=1;
                        while n<=Ngoods && ~isempty(locations)
                            [Y,I]=sort(g.L(locations),'descend'); % sort by population
                            param.Zjn(locations(I(1)),n+1)=1; % only largest cities produce one of the specific goods
                            param.Zjn(locations(I(1)),1)=0;
                            
                            % remove largest city and its neighbors from the list
                            locations=setdiff(locations,[locations(I(1)),g.nodes{locations(I(1))}.neighbors]);
                            n=n+1;
                        end
                        param.N=1+n-1;
                        param.Zjn = param.Zjn(:,1:param.N);
                        param.m = param.m(1:param.N);
                        
                        %identify differentiated nodes
                        nodes_dif = find(sum(param.Zjn(:,2:end)>0,2));
                        
                    case 'nuts'
                        
                        % This version assigns one good to largest city in
                        % each nuts for the Ngoods largest nuts only
                        
                        % number of differentiated goods
                        Ngoods = min( nb_NUTS,Ngoods_fixed );
                        
                        % largest NUTS
                        pop_NUTS = zeros( 1,nb_NUTS );
                        for i=1:nb_NUTS
                            pop_NUTS(i) =sum( g.L( g.region==NUTS_list(i) ) ); % total population of nuts
                        end
                        [ ~,diff_NUTS_index ] = sort(pop_NUTS,'descend');
                        
                        % productivity matrix
                        param.Zjn = zeros( g.J,Ngoods+1 );
                        param.Zjn(:,1) = 1; % all cities produce the agricultural good
                        
                        n=1;
                        for i=diff_NUTS_index( 1:Ngoods )          % the differentiated goods are ordered based on the ranking of the nuts that produces it
                            
                            list_nodes = find( g.region==NUTS_list(i) );
                            [~,largest_id]=max( g.L(list_nodes) ); % identify most populated cell in region
                            
                            param.Zjn( list_nodes( largest_id ),n+1 )=1; % only largest cell in region produces differentiated goods
                            param.Zjn( list_nodes( largest_id ),1 )=0;   % and doesn't produce homogeneous good
                            n=n+1;
                            
                        end
                        
                        param.N=1+Ngoods;
                        param.m = ones(param.N,1);
                        
                        %identify differentiated nodes
                        nodes_dif = find(sum(param.Zjn(:,2:end)>0,2));
                        
                end
                
                %% calibrate
                verbose=true;
                calibrate();
                results_actual = results; % keep the calibrated allocation for welfare comparison
                results_actual.Ijk = g.avI;
                if calibrate_d0==true && save_calibrated_d0==true
                    save_d0_to_file( d0_file,beta,gamma,param.mobility,param.cong,Ngoods_fixed,d0 );
                end
                
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
                
                %% project delta_i_tilde on fundamentals
                
                % construct dataset at bilateral level
                id=find( tril(g.adjacency)==1 );
                Npairs = length(unique_edges);  % number of unique (unordered) origin-destination pairs
                
                d_i_tilde = zeros( Npairs,1 );
                dist = zeros( Npairs,1 );
                diff_alt = zeros( Npairs,1 );
                av_alt = zeros( Npairs,1 );
                av_rugged = zeros( Npairs,1 );
                %FE = zeros( Npairs,g.J );  % fixed effect
                
                data = cell( Npairs,10 );
                
                for i=1:Npairs   % each row is a pair
                    
                    node1 = unique_edges( i,1 );
                    node2 = unique_edges( i,2 );
                    
                    d_i_tilde( i ) = log( delta_i_tilde( node1,node2 ) );
                    dist( i ) = log( g.distance( node1,node2 ) );
                    diff_alt( i ) = abs( log( g.av_alt( node1 )/g.av_alt( node2 ) ) );
                    av_alt( i ) = 1/2*( log( g.av_alt( node1 ) )+log( g.av_alt( node2 ) ) );
                    av_rugged( i ) = 1/2*( log( g.rugged( node1 ) )+log( g.rugged( node2 ) ) );
                    
                    % matrix to export to csv
                    data{i,1} = country;
                    data{i,2} = country_icc;
                    data{i,3} = node1;
                    data{i,4} = node2;
                    data{i,5} = delta_i_tilde( node1,node2 );
                    data{i,6} = g.distance( node1,node2 );
                    data{i,7} = g.av_alt( node1 );
                    data{i,8} = g.av_alt( node2 );
                    data{i,9} = g.rugged( node1 );
                    data{i,10} = g.rugged( node2 );
                    
                end
                
                write_to_excel = 0;
                if write_to_excel
                    
                    % export to excel
                    headers = { 'Country','ICC','Cell1','Cell2','dtilde','dist','altitude1','altitude2','rugged1','rugged2' };
                    xlswrite( [ path_save_tables,filename,'_dtilde.xls' ],headers,country_icc,'A1' );
                    xlswrite( [ path_save_tables,filename,'_dtilde.xls' ],data,country_icc,'A2' );
                    
                end
                
                % run linear model
                tbl = table( d_i_tilde,dist,diff_alt,av_alt,av_rugged,...
                    'VariableNames',{'DeltaTilde','Distance','Diff_Alt','Av_Alt','Av_Rugged'} );
                mdl = fitlm( tbl,'DeltaTilde~Distance+Diff_Alt+Av_Rugged' );
                
                % project delta_i
                g.delta_i_projected = zeros(g.J,g.J);
                for k=1:Npairs
                    g.delta_i_projected(unique_edges(k,1),unique_edges(k,2)) = exp(mdl.Fitted(k));
                end
                g.delta_i_projected=g.delta_i_projected+g.delta_i_projected';
                g.delta_i_projected( g.adjacency == 0 ) = 0;
                g.delta_i_projected = g.delta_i_projected / sum( g.delta_i_projected(:).*g.avI(:) );  % normalize to meet constraint with K=1
                
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
end
toc

