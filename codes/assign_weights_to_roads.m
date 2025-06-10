% This code: assign a measure of road infrastructure to actual road network data

function roads = assign_weights_to_roads( roads )

% keep roads that are not under construction
if isfield(roads,'EXS')
    roads = roads( cell2mat({roads.EXS})==28 );
end

% number of lanes
road_lanes_mat = cell2mat( {roads.LTN} );       % number of lanes
road_lanes_missing = ( road_lanes_mat<=0 );     % missing lanes
road_lanes_mat( road_lanes_mat<=0 ) = 1;        % replace missing number of lanes with 1 lane

% max_road_lanes=max( road_lanes_mat );           % max number of road lanes
% tabulate(road_lanes_mat)

% intended use  % =0: local; =1: national
% roads = roads( cell2mat({roads.RTT})~=0 ); % keep roads with non-missing intended use
road_use_mat = cell2mat( {roads.RTT} );
road_use_missing = ( road_use_mat<=0 );     % missing road use
road_use_mat( road_use_mat<=0 ) = 0;      % unknown
road_use_mat( road_use_mat==15 ) = 0;     % secondary
road_use_mat( road_use_mat==984 ) = 0;    % local
road_use_mat( road_use_mat==14 ) = 0;     % primary 
road_use_mat( road_use_mat==16 ) = 1;     % national
%tabulate( road_use_mat )

% median?   % =0: without median =1: with median; 
road_median_mat = cell2mat( {roads.MED} );
road_median_mat(road_median_mat<=0 ) = 2; % unknown-->no median
road_median_mat(road_median_mat==2 ) = 0; % no median = 0
%tabulate(road_median_mat);

% paved?  % =0: unpaved; =1: paved
road_paved_mat = cell2mat( {roads.RST} );
road_paved_mat(road_paved_mat<=0 ) = 2;  % unknown-->unpaved
road_paved_mat(road_paved_mat==2 ) = 0;
%tabulate(road_paved_mat);

% distance
road_length_mat = deg2km( cell2mat( {roads.Shape_Leng} ) );

% cheapest-path coefficients
chi_lane = 0.2;
chi_use = 1.07;
chi_paved = 1.35;
chi_median = 1.05;

% weight
shape_weight = road_length_mat.*...
                ( road_lanes_mat.^( -chi_lane ) ).*...
                 chi_use.^( 1-road_use_mat ).*...
                 chi_paved.^( 1-road_paved_mat ).*...
                 chi_median.^( 1-road_median_mat );   % weight of entire shape

%% assign a distance in km and a weight to each edge

for j=1:length(roads)
    
    [ j length(roads) ]
    
    roads(j).totdist = road_length_mat(j); % deg2km( roads(j).SHAPE_Leng ); % total distance, in KM
    if isempty(roads(j).totdist)
        roads(j).totdist=0; 
    end
    
    % length of each link in km
    Ntemp = length(roads(j).X)-1;
    x0 = roads(j).X(1:Ntemp-1);
    y0 = roads(j).Y(1:Ntemp-1);
    x1 = roads(j).X(2:Ntemp);
    y1 = roads(j).Y(2:Ntemp);  
    if roads(j).totdist>0
        roads(j).distances = distance( y0,x0,y1,x1 )./sum( distance( y0,x0,y1,x1 ) ) * roads(j).totdist;       
    elseif roads(j).totdist==0
        roads(j).distances = distance( y0,x0,y1,x1 );
    end
    
    % link features
    roads(j).lanes = road_lanes_mat(j);
    roads(j).totlanes = roads(j).lanes*roads(j).totdist;    
    roads(j).use = road_use_mat(j);
    roads(j).median = road_median_mat(j);
    roads(j).paved = road_paved_mat(j);
    
    % missing values
    roads(j).missing_lanes = road_lanes_missing(j);
    roads(j).missing_use = road_use_missing(j);
    
    % segment weight
    roads(j).weights = distance( y0,x0,y1,x1 )./sum( distance( y0,x0,y1,x1 ) ) * shape_weight(j);
    
end