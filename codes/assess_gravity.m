%%% assess gravity properties for Spain

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

% parameters of grid to load
x_ver_hor = 0.6;   
x_diag = 0.6;      

% parameters that remain constant throughout the loop across calibrations
alpha = 0.4;
beta = 0.13;
gamma = 0.10;
sigma = 5;
rho = 0;
a = 1;
nu = 1;
mobil = 1;
cong = 1;

% number of differentiated goods
Ngoods = 15;  % for 'largest ciies', this is the number of goods
              % for 'nuts', the number of goods is the min between Ngoods and the number of NUTS regions
              % this may change if we use 'nuts' and the number of nuts is less than Ngoods

% how to allocate cities?
city_allocation = 'nuts'; % 'largest cities' or 'nuts'              
NUTS = 2;

%% load Spain
country = 'Spain';

% calibration
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
                            '_mobil',num2str( mobil ),...
                            '_cong',num2str( cong ),...
                            '_ngoods',num2str( Ngoods ),...
                            '_EGM',num2str(EGM) ];

if strcmp(city_allocation,'nuts')                 
    
    filename = [ filename,'_nuts',num2str(NUTS) ];

end                                              
                        
load( [ path_save_calibrations,filename,'_calib.mat' ] ) % load calibration structure

% Unpack
g = calibration.g;
param = calibration.param;
results = calibration.results_actual;

% country graph
graph_filename = [ path_load_grids,country,'_grid_',...
                               num2str( x_diag ),'_',...
                               num2str( x_ver_hor ),'_EGM8.mat' ];
load( [ graph_filename ] );                       

g.region = [ country_graph.graph_export.nuts.codes{:,NUTS} ]';                                              
bilateral_dist = country_graph.graph_export.nuts_dist_matrix{NUTS};

%% generate trade flows from calibrated model

[ intra_reg_trade_share_differentiated,gravity_coeff,bilateral_trade_mat ] = compute_trade_shares_and_gravity(param,g,results,1,bilateral_dist);

%% clean up names

% mapping from region to region id
    
% region names
[~,I,~] = unique( cell2mat(g.nodes_nuts(:,1)) );
nuts2_list(:,1) = g.nodes_nuts(I,1);
nuts2_list(:,2) = g.nodes_nuts(I,3);

names_clean{1}='Galicia';
names_clean{2}='Asturias';
names_clean{3}='Cantabria';
names_clean{4}='Pais Vasco';
names_clean{5}='Navarra';
names_clean{6}='La Rioja';
names_clean{7}='Aragon';
names_clean{8}='Madrid';
names_clean{9}='Castilla Leon';
names_clean{10}='Castilla La Mancha';
names_clean{11}='Extremadura';
names_clean{12}='Cataluna';
names_clean{13}='Valencia';
names_clean{14}='Andalucia';
names_clean{15}='Murcia';

% generate long dataset
bilateral_trade_calib = bilateral_trade_mat;
bilateral_trade_data = country_graph.trade_matrix.nuts2;
bilateral_dist = country_graph.dist_matrix.nuts2;

trade_dataset_long = {};
j=1;
nb_NUTS = size(bilateral_dist,1);
for o=1:nb_NUTS
    for d=1:nb_NUTS
        var_n = 1;  % order of variable in dataset
        trade_dataset_long{j,var_n}=o; var_n=var_n+1; % use original NUTS to construct the trade dataset
        trade_dataset_long{j,var_n}=names_clean{o}; var_n=var_n+1;
        trade_dataset_long{j,var_n}=d; var_n=var_n+1;
        trade_dataset_long{j,var_n}=names_clean{d}; var_n=var_n+1;
        trade_dataset_long{j,var_n}=bilateral_dist(o,d); var_n=var_n+1;
        trade_dataset_long{j,var_n}=bilateral_trade_data(o,d); var_n=var_n+1;
        trade_dataset_long{j,var_n}=bilateral_trade_calib(o,d);
        j=j+1;
    end
end
%% compute gravity

write_excel = 1;
if write_excel
    headers = {'origin','origin_name',...
               'destination','destination_name',...
               'distance','x_real','x_calib'};
    switch city_allocation
        case 'largest cities'
            xlswrite( [ path_gravity,'gravity_Spain_largest_cities_Ngoods',num2str(Ngoods),'.xls' ],headers,'trade_Spain','A1' );
            xlswrite( [ path_gravity,'gravity_Spain_largest_cities_Ngoods',num2str(Ngoods),'.xls' ],trade_dataset_long,'trade_Spain','A2' );
        case 'nuts'
            xlswrite( [ path_gravity,'gravity_Spain_nuts_Ngoods',num2str(Ngoods),'.xls' ],headers,'trade_Spain','A1' );
            xlswrite( [ path_gravity,'gravity_Spain_nuts_Ngoods',num2str(Ngoods),'.xls' ],trade_dataset_long,'trade_Spain','A2' );
    end
end                      