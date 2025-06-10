%% save coordinates corresponding to country boundaries

clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;

%% load country list

country_list_short;

%% create polygon with country boundaries for each country and their union

Xmin_all_countries = Inf;
Xmax_all_countries = -Inf;
Ymin_all_countries = Inf;
Ymax_all_countries = -Inf;

all_countries_bounds = [];

%%
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
    
     if strcmp(country,'Denmark') || strcmp(country,'Germany') || strcmp(country,'Spain') || ...
           strcmp(country,'Finland') || strcmp(country,'Bulgaria') || strcmp(country,'Croatia') || ...
           strcmp(country,'Greece')  || strcmp(country,'Netherlands') || strcmp(country,'Norway') ...
           || strcmp(country,'Sweden')

            %international boundaries
            file = [datafolder_ne,'ne_10m_admin_0_map_units.shp'];
            country_bounds_temp = shaperead(file,'Selector',{@(name)strcmp( name,country ),'GEOUNIT'}); 
            
            X = country_bounds_temp.X;
            Y = country_bounds_temp.Y;

     elseif strcmp(country,'Serbia')

                    if ~ispc()
                        datafolder_serbia = [datafolder_DIVA,'/SRB_adm/'];
                    else
                        datafolder_serbia = [datafolder_DIVA,'\SRB_adm\'];
                    end

                    %international boundaries
                    country_bounds_temp = shaperead([datafolder_serbia,'SRB_adm0.shp']);

                    X = country_bounds_temp.X;
                    Y = country_bounds_temp.Y; 

     else   
    
        %international boundaries
        file  = [datafolder_Eurogeog_country,'PolbndL.shp'];     
        country_bounds_temp = shaperead( file,'Selector',{@(use)use==23,'USE'} );  % international boundaries
    
        % in EGM, country bounds are lines, define country bound as a polygon   
        [ X,Y ] = polymerge( cell2mat({country_bounds_temp.X}),cell2mat({country_bounds_temp.Y}) );     
        [ X,Y ] = poly2cw( X,Y );   

     end           

     % set country bounds                                  
        country_bounds.Geometry = 'Polygon';
        country_bounds.BoundingBox = [ min( X( ~isnan(X) ) ) min( Y( ~isnan(Y) ) );...
                                       max( X( ~isnan(X) ) ) max( Y( ~isnan(Y) ) ) ];                                 
        country_bounds.X = X;
        country_bounds.Y = Y;        
     
     
     %% Chosen region may include islands. Keep largest polygon
     
     [lat,lon] = polysplit( country_bounds.Y',country_bounds.X' ); % split the country into its polygons
     
     if length(lat)>1
         polyareas = [];
         for j=1:length(lat)
             polyareas(j) = polyarea( lon{j},lat{j} );
        end

         if strcmp(country,'Denmark')
            [~,keep_poly] = maxk( polyareas,3 ); % for Denmark, keep largest islands too
         else
            [~,keep_poly] = maxk( polyareas,1 );
         end

         [ country_bounds.Y,country_bounds.X ] = polyjoin( lat(keep_poly),lon(keep_poly) );
     end
     
     %% replace bounds with largest polygons
     X = cell2mat({country_bounds.X});
     Y = cell2mat({country_bounds.Y});
     country_bounds.BoundingBox = [ min( X( ~isnan(X) ) ) min( Y( ~isnan(Y) ) );...
                                    max( X( ~isnan(X) ) ) max( Y( ~isnan(Y) ) ) ];  
     % add name
     country_bounds.country = country;
     
     % save country_bounds
     file = [ datafolder_Eurogeog_country,'country_bounds.mat'];
     save(file,'country_bounds');
     
     % record the min and max coordinates across all countries
     Xmin_all_countries = min( Xmin_all_countries,country_bounds.BoundingBox(1,1) );
     Ymin_all_countries = min( Ymin_all_countries,country_bounds.BoundingBox(1,2) );
     Xmax_all_countries = max( Xmax_all_countries,country_bounds.BoundingBox(2,1) );
     Ymax_all_countries = max( Ymax_all_countries,country_bounds.BoundingBox(2,2) );
     
     BoundingBox_allcountries = [ Xmin_all_countries Ymin_all_countries;...
                                  Xmax_all_countries Ymax_all_countries ];                             
     
     % add to file with all countries
     all_countries_bounds(country_n).Geometry = 'Polygon';
     all_countries_bounds(country_n).X = country_bounds.X;
     all_countries_bounds(country_n).Y = country_bounds.Y;
     all_countries_bounds(country_n).country = country_bounds.country;
     all_countries_bounds(country_n).BoundingBox = country_bounds.BoundingBox;
     
    %% plot
    h=figure;
    mapshow(country_bounds,'FaceColor','none')
    title(country) 
    saveas(h, [path_save_grids,'country_bounds/',country,'_bounds.jpg']) 
     
end

%% save

% save file for all countries
file = [ datafolder_Eurogeog,'all_countries_bounding_box.mat'];
save(file,'BoundingBox_allcountries');

file = [ datafolder_Eurogeog,'all_countries_bound.mat'];
save(file,'all_countries_bounds');

%% save the boundaries of the connected countries

clear
close all
clc

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;

% load the set of connected countries
file = [ datafolder_Eurogeog,'all_countries_bound.mat'];
load(file);

% define a region that is the union of connected countries
drop = {'Cyprus','Finland','Georgia','Ireland','Lithuania','Latvia','Macedonia','Northern Ireland','Moldova',...
        'Slovenia','Hungary','Slovakia','Czech Republic'};
%drop = {'Cyprus','Finland','Georgia','Ireland','Lithuania','Latvia','Macedonia','Northern Ireland','Moldova'};

    
for i=1:length(all_countries_bounds)
        keep(i) = 1-max( strcmp( all_countries_bounds(i).country,drop ) );
end

h=figure;
country_bounds = all_countries_bounds( logical(keep) );
mapshow( country_bounds,'FaceColor','none' )
title('Connected Countries')

% save
file = [ datafolder_connected_countries,'connected_countries_bound_many.mat'];
save(file,'country_bounds');  % this structure contains one entry for each country

% plot
h=figure;
mapshow(all_countries_bounds,'FaceColor','none')
mapshow(country_bounds,'FaceColor','red')
title('All countries and connected countries')  
saveas(h, [path_save_grids_country_bounds,'all_countries_bounds.jpg']) 

% set the union of all countries
X = country_bounds(1).X;
Y = country_bounds(1).Y;
for i=2:length(country_bounds)
    [X,Y] = polybool('union',X,Y,...
                      country_bounds(i).X,country_bounds(i).Y);
end

clear country_bounds
country_bounds.Geometry = 'Polygon';
country_bounds.X = X;
country_bounds.Y = Y;
country_bounds.BoundingBox = [ min( X( ~isnan(X) ) ) min( Y( ~isnan(Y) ) );...
                                       max( X( ~isnan(X) ) ) max( Y( ~isnan(Y) ) ) ];
% save
file = [ datafolder_connected_countries,'connected_countries_bound_single.mat'];
save(file,'country_bounds');  % this structure contains a single entry with all countries