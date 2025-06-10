%% split ETOPO into countries and into connected countries in Europe
% compute ruggedness measure within 0.25 degree cells

clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')
warning('off')

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;

%% load country list

country_list_short;

%% load Europe grid from ETOPO

% load TIFF
[A,R]  = geotiffread( [datafolder_ETOPO,....
                        'ETOPO1_Ice_g_geotiff'] ); 

%{
  This is a GeographicCellsReference with properties:

             LatitudeLimits: [1x2 double]
            LongitudeLimits: [-180 180.000000000002]
                 RasterSize: [17400 43200]
       RasterInterpretation: 'cells'
           ColumnsStartFrom: 'north'
              RowsStartFrom: 'west'
       CellExtentInLatitude: 0.00833333333333339
      CellExtentInLongitude: 0.00833333333333339
     RasterExtentInLatitude: 145.000000000001
    RasterExtentInLongitude: 360.000000000002
           XIntrinsicLimits: [0.5 43200.5]
           YIntrinsicLimits: [0.5 17400.5]
       CoordinateSystemType: 'geographic'
                  AngleUnit: 'degree'
%}

% Locate NW corners of each cell

    % Y runs North-South
    Ygrid = R.YWorldLimits(2):-R.CellExtentInWorldY:R.YWorldLimits(1)+R.CellExtentInWorldY;

    % X runs West-East
    Xgrid = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2)-R.CellExtentInWorldX;

% keep cells in a box containing all the countries in our data
    file = [ datafolder_Eurogeog,'all_countries_bounding_box.mat'];
    load(file);
    minX = floor( BoundingBox_allcountries(1,1) );
    maxX = ceil( BoundingBox_allcountries(2,1) );
    minY = floor( BoundingBox_allcountries(1,2) );
    maxY = ceil( BoundingBox_allcountries(2,2) );    

    keepX = logical( ( Xgrid>=minX ).*( Xgrid<=maxX ) );
    keepY = logical( ( Ygrid>=minY ).*( Ygrid<=maxY ) );

    Ygrid = Ygrid( keepY );
    Xgrid = Xgrid( keepX );
    A = A( keepY,keepX );
    
% save
    etopodata.A = A;
    etopodata.Xgrid = Xgrid;
    etopodata.Ygrid = Ygrid;
    etopodata.R = R;

    save([datafolder_ETOPO,'europe_etopo.mat'],'etopodata');

%% construct measures of ruggedness, aggregate into 0.25 degree cells and generate map structure

% set cell size, in minutes
cellsize = 0.25;
NN = round( cellsize/R.CellExtentInWorldY );  % each 0.5 cell will feature NN^2 cells from the map

% north-west corners of each cell
Ygrid_coarse = round(max(Ygrid)):-cellsize:min(Ygrid);
Xgrid_coarse = round(min(Xgrid)):cellsize:max(Xgrid);

clc
counter = 1;

for j=1:length(Xgrid_coarse)
    
    [ j length(Xgrid_coarse) ]
    
    for i=1:length(Ygrid_coarse)
        
        places(counter).Geometry = 'Polygon';
        places(counter).X = [ Xgrid_coarse(j)          Xgrid_coarse(j) Xgrid_coarse(j)+cellsize Xgrid_coarse(j)+cellsize NaN ];
        places(counter).Y = [ Ygrid_coarse(i)-cellsize Ygrid_coarse(i) Ygrid_coarse(i)          Ygrid_coarse(i)-cellsize NaN ];
        
        places(counter).centerX = mean( places(counter).X( ~isnan( places(counter).X ) ) );
        places(counter).centerY = mean( places(counter).Y( ~isnan( places(counter).Y ) ) );
        
        temp = A( NN*( i-1 )+1:NN*i,...
                  NN*( j-1 )+1:NN*j );
              
        % compute differences in elevation
        x = temp(:,2:NN)-temp(:,1:NN-1);
        diff_elev = reshape( x,size( x,1 )*size( x,2 ),1 ); %horizontal
        
        x = temp(2:NN,:)-temp(1:NN-1,:); 
        diff_elev = [ diff_elev;reshape( x,size( x,1 )*size( x,2 ),1 ) ];  %vertical
        
        for row = 0:NN-2
            x = diag( temp( 1:NN-row,1+row:NN ) );     
            diff_elev = [ diff_elev;x ];
        end
        for row = 1:NN-2
            x = diag( temp( 1+row:NN,1:NN-row ) );     
            diff_elev = [ diff_elev;x ];
        end
        
        temp = temp(:,NN:-1:1); % reorder columns
        for row = 0:NN-2
            x = diag( temp( 1:NN-row,1+row:NN ) );     
            diff_elev = [ diff_elev;x ];
        end
        for row = 1:NN-2
            x = diag( temp( 1+row:NN,1:NN-row ) );   
            diff_elev = [ diff_elev;x ];
        end
              
        places(counter).av_alt = mean( temp(:) );             % average altitude
        places(counter).sd_alt = std( double( temp(:) ) );    % sd of altitude
        places(counter).rugged = std( double( diff_elev ) );  % ruggedness index
        
        counter = counter+1;
        
    end
end

                    
%% retain the section of European grid corresponding to each country
for country_n=1:Ncountries
    
    country_icc = char( countries(country_n) );  % country ICC code
    country = icc2name( country_icc );           % country name

    % load country_bounds
    if ~ispc()   %% we are in a mac

        % folder to save the outcome
        datafolder_Eurogeog_country =  [ datafolder_Eurogeog,country_icc,'/'];

    else   %% we are in a pc   

        % folder to save the outcome
        datafolder_Eurogeog_country =  [ datafolder_Eurogeog,country_icc,'\' ];

    end
    
    clc
    disp( ['Country: ',country] )     
    
    % load bounds 
    file  = [ datafolder_Eurogeog_country,'country_bounds.mat'];
    load(file);
    
    % keep places within the defined boundaries
    keep_place = inpolygon( cell2mat({places.centerX}'), cell2mat({places.centerY}'),...
                          [ min(country_bounds.X)-1 max(country_bounds.X)+1 ],...
                          [ min(country_bounds.Y)-1 max(country_bounds.Y)+1 ] );
    relief_etopo = places( keep_place );
    
    % save
    file  = [ datafolder_Eurogeog_country,'relief_etopo.mat'];
    save(file,'relief_etopo');

    % plot
    %{
    h=figure;
    mapshow(relief_etopo,'FaceColor','white')
    mapshow(country_bounds,'FaceColor','none')
    title(country)
    if ismac
        saveas(h, [datafolder_ETOPO,'maps/',country,'_map_etopo.pdf']) 
    else
        saveas(h, [datafolder_ETOPO,'maps\',country,'_map_etopo.pdf']) 
    end
    %}
       
end

%% combine into one file for the connected countries in Europe

clear
close all
clc

% load list of connected countries
folders;
country_list_connected;

% append the files
for country_n=1:Ncountries
    
    country_icc = char( countries_connected (country_n) );  % country ICC code
    country = icc2name( country_icc );           % country name

    if ~ispc()   %% we are in a mac

        % folder to save the outcome
        datafolder_Eurogeog_country = [ datafolder_Eurogeog,country_icc,'/'];

    else   %% we are in a pc   

        % folder to save the outcome
        datafolder_Eurogeog_country = [ datafolder_Eurogeog,country_icc,'\' ];

    end
    
    clc
    disp( ['Country: ',country] )    
    
    % load sedac
    file  = [ datafolder_Eurogeog_country,'relief_etopo.mat'];
    temp = load(file); 
    
    % save temp file
    relief_etopo_temp{country_n} = temp.relief_etopo;
    for i=1:length(relief_etopo_temp{country_n})
        relief_etopo_temp{country_n}(i).country = country;
    end
    
end

%
i=1;
for country_n=1:Ncountries
    for j=1:length(  relief_etopo_temp{country_n} )
        relief_etopo(i).Geometry =  relief_etopo_temp{country_n}(j).Geometry;
        relief_etopo(i).X =  relief_etopo_temp{country_n}(j).X;
        relief_etopo(i).Y =  relief_etopo_temp{country_n}(j).Y;
        relief_etopo(i).centerX =  relief_etopo_temp{country_n}(j).centerX;
        relief_etopo(i).centerY =  relief_etopo_temp{country_n}(j).centerY;
        relief_etopo(i).av_alt =  relief_etopo_temp{country_n}(j).av_alt;
        relief_etopo(i).sd_alt =  relief_etopo_temp{country_n}(j).sd_alt;
        relief_etopo(i).rugged =  relief_etopo_temp{country_n}(j).rugged;
        relief_etopo(i).country =  relief_etopo_temp{country_n}(j).country;
        i=i+1;
    end
end

% save
save([datafolder_connected_countries,'relief_etopo.mat'],'relief_etopo');

% load connected country bounds
file = [ datafolder_connected_countries,'connected_countries_bound_many.mat'];
load(file);

% plot    
%{
h=figure
Colors = makesymbolspec('Polygon', {'av_alt', ...
    [0 max( [relief_etopo.av_alt] )],...
    'FaceColor', jet(numel(relief_etopo)) });
mapshow(relief_etopo,'SymbolSpec',Colors) 
mapshow(country_bounds,'FaceColor','none','EdgeColor','white')  
saveas(h, [datafolder_connected_countries,'map_etopo.pdf']);
%}