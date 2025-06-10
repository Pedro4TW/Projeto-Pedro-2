clear
close all
clc

%% preliminaries on the raw data, country by country and for connected countries in Europe

% these steps create country_graph.m called by the calibration

    % determine country bounds
    define_country_boundaries;

    % aggregate the population data to small cells within each country
    split_SEDAC;

    % break down income data into individual countries
    split_gecon;

    % break down topography data into individual countries
    split_ETOPO;

    % generate country grids and actual discretized networks
    make_country_graphs;
    make_country_graphs_largecells; % same as above, but all cells are size 1 by 1 
    make_country_graphs_smallcells; % same as above, but all cells are size 0.5 by 0.5

    % make_connected_countries graph
    make_connected_countries_graphs; 
    
    % allocate nodes to NUTS2 regions for all countrie available NUTS classification
    allocate_nodes_to_pol_regions_all;  % run three times: for larce_cells=1 and =0, and for small_cells=1
    
    % allocate Spanish cells to Spanish regions, used to assess gravity
    allocate_nodes_to_pol_regions_spain;

%% calibration
    run_calibration; % calibrate all 24 countries separately       
    calibrate_europe; % calibrate the european case

%% counterfactuals
    run_counterfactuals;
    run_counterfactuals_europe;

%% prepare outcomes for tables
    make_summary_tables_revision;
    assess_gravity;

%% figures and tables
         
% make_figures_2_3_5.m     % figures 2,3,5
% make_figure_6            % figures 6
% make_figure_A2           % figure A2
% figure4.do               % figure 4
% make_table_A1.m          % tableA1
% all_tables.do            % remaining tables and figures