%% split SEDAC into countries and into connected countries in Europe
% aggregate the population data to evenly spaced points at the center of 0.1 degree cells
% each point is assigned the population of its cell
% population is further aggregated to 0.5 degree cells by make_country_graps.m

clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;

%% load country list

country_list_short;

%% load Europe grid from SEDAC

% load TIF, this is a 250m cell grid
[ A,R ]  = geotiffread( [datafolder_SEDAC,....
                        'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2005.tif'] );
                   
% zeros appear as tiny negative numbers but they add up, set them to zero                    
A( A<0 )=0;

% Locate NW corners of each cell

    % Y runs North-South
    Ygrid = R.LatitudeLimits(2):-R.CellExtentInLatitude:R.LatitudeLimits(1)+R.CellExtentInLatitude;

    % X runs West-East
    Xgrid = R.LongitudeLimits(1):R.CellExtentInLongitude:R.LongitudeLimits(2)-R.CellExtentInLongitude;

% keep cells in a box containing all countries in the data
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
    sedac.A = A;
    sedac.Xgrid = Xgrid;
    sedac.Ygrid = Ygrid;
    sedac.R = R;

    save([datafolder_SEDAC,'europe_sedac.mat'],'sedac');

%% aggregate into coarser cells and generate places map structure

% set cell size, in minutes
cellsize = 0.1;
NN = round( cellsize/R.CellExtentInLatitude ); % each cell in coarser grid will contain NN^2 cells in finer original grid

% centers of cells
Ygrid_coarse = max(Ygrid)-cellsize/2:-cellsize:min(Ygrid);
Xgrid_coarse = min(Xgrid)+cellsize/2:cellsize:max(Xgrid);

clc
counter = 1;
for j=1:length(Xgrid_coarse)
    
    [ j length(Xgrid_coarse) ]
    
    for i=1:length(Ygrid_coarse)
        
        places(counter).Geometry = 'Point';
        places(counter).X = Xgrid_coarse(j);
        places(counter).Y = Ygrid_coarse(i);
        
        temp = A( NN*( i-1 )+1:NN*i,...
                  NN*( j-1 )+1:NN*j );
              
        places(counter).population = sum( temp(:) );  % the center of each cell is assigned its population
        
        counter = counter+1;
        
    end
end
                    
%% retain the section of grid corresponding to each country

for country_n=1:Ncountries
    
    country_icc = char( countries(country_n) );  % country ICC code
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
    
    % load country_bounds
    file  = [ datafolder_Eurogeog_country,'country_bounds.mat']; 
    load(file);
    
    % keep places within the defined boundaries
    keep_pop = inpolygon( cell2mat({places.X}'), cell2mat({places.Y}'), country_bounds.X, country_bounds.Y );
    places_sedac = places( keep_pop );
    
    % save place
    file  = [ datafolder_Eurogeog_country,'popplaces_sedac.mat'];
    save(file,'places_sedac');
    
    %% plot
    %{
    h=figure;
    mapshow(places_sedac,'Marker','.','MarkerSize',10);
    mapshow(country_bounds,'FaceColor','none')
    totpop =sum( cell2mat( {places_sedac.population} ) ) ;
    title([country,'. Tot Pop (Million): ',num2str(round(totpop/10000)/100)]); 
    if ismac
        saveas(h, [datafolder_SEDAC,'maps/',country,'_map_sedac.pdf']) 
    else
        saveas(h, [datafolder_SEDAC,'maps\',country,'_map_sedac.pdf']) 
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
    file  = [ datafolder_Eurogeog_country,'popplaces_sedac.mat'];
    temp = load(file); 
    
    % save temp file
    places_sedac_temp{country_n} = temp.places_sedac;
    for i=1:length(places_sedac_temp{country_n})
        places_sedac_temp{country_n}(i).country = country;
    end
    
end

i=1;
for country_n=1:Ncountries
    for j=1:length( places_sedac_temp{country_n} )
        places_sedac(i).Geometry =  places_sedac_temp{country_n}(j).Geometry;
        places_sedac(i).X =  places_sedac_temp{country_n}(j).X;
        places_sedac(i).Y =  places_sedac_temp{country_n}(j).Y;
        places_sedac(i).population =  places_sedac_temp{country_n}(j).population;
        places_sedac(i).country =  places_sedac_temp{country_n}(j).country;
        i=i+1;
    end
end

% save
save([datafolder_connected_countries,'popplaces_sedac.mat'],'places_sedac');

% load connected country bounds
file = [ datafolder_connected_countries,'connected_countries_bound_many.mat'];
load(file);

%% plot    
%{
h=figure;
mapshow(places_sedac,'Marker','.','MarkerSize',10);
mapshow(country_bounds,'FaceColor','none')
totpop =sum( cell2mat( {places_sedac.population} ) );
saveas(h, [datafolder_connected_countries,'connected_countries.jpg']);
%}