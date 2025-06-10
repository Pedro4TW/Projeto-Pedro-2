%% split data from G-econ into countries and into connected countries in Europe
% identify each cell in Gecon with a country in our data

clear
close all
clc

set(0,'DefaultFigureWindowStyle','docked','DefaultFigureVisible','on')

% -----------------------
% DEFINE FOLDER LOCATIONS
folders;

%% transform G-Econ dataset into map structure
gecon_data = csvread([datafolder_gecon,'gecon_data.csv']);

lat = gecon_data(:,1);
lon = gecon_data(:,2);
pop = gecon_data(:,3);
Y = gecon_data(:,4);

N = length(Y);

countries = importdata([datafolder_gecon,'gecon_countries.csv']);

for j=1:N
    
    places_gecon(j).Geometry = 'Polygon';
   
    places_gecon(j).X = [ lon(j) lon(j)   lon(j)+1 lon(j)+1  NaN];   % lon(j), lat(j) is the SW corner
    places_gecon(j).Y = [ lat(j) lat(j)+1 lat(j)+1 lat(j)    NaN];

    places_gecon(j).Xcorner = lon(j); % lon(j), lat(j) is the SW corner
    places_gecon(j).Ycorner = lat(j);
    
    places_gecon(j).population = pop(j);
    places_gecon(j).gdp = Y(j);
    places_gecon(j).country = countries(j);
    
end

save([datafolder_gecon,'places_gecon'])

%% load country list

country_list_short;

%% loop through all countries

for country_n=1:Ncountries

    %%
    country_icc = char( countries(country_n) );  % country ICC code
    country = icc2name( country_icc );           % country name

    if ~ispc()   %% we are in a mac

        % folder to save the outcome
        datafolder_Eurogeog_country =  [ datafolder_Eurogeog,country_icc,'/'];

    else   %% we are in a pc   

        % folder to save the outcome
        datafolder_Eurogeog_country =  [ datafolder_Eurogeog,country_icc,'\' ];

    end
    
    clc
    disp( ['Country: ',country] )   
    
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% bring income data %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % bring all G-Econ
    file  = [datafolder_gecon,'places_gecon.mat'];
    temp = load(file);
    places_gecon = temp.places_gecon;
      
    % keep a cell in G-Econ if:
    % it has the same country name as the chosen country
    
    keep = zeros( length(places_gecon),1 );
    if ~strcmp( country_icc,'ND')
        for j=1:length(places_gecon)
            keep(j) = strcmp( places_gecon(j).country,country );
        end
        places_gecon = places_gecon( logical(keep) );
    elseif strcmp(country_icc,'ND')   % Sections in Northern Ireland appears as Great Britain 
        for j=1:length(places_gecon)
            keep(j) = strcmp( places_gecon(j).country,'Great Britain' );
        end
        places_gecon = places_gecon( logical(keep) );
    end
    
    % further restrict to cells that overlap with country boundaires
    file  = [ datafolder_Eurogeog_country,'country_bounds.mat'];
    load(file);
    keep = zeros( length(places_gecon),1 );
    for j=1:length(places_gecon)
          keep(j) = ~isempty( polybool( 'intersection',...
                              cell2mat({places_gecon(j).X}), cell2mat({places_gecon(j).Y}),...
                              country_bounds.X, country_bounds.Y ) );

     end
     places_gecon = places_gecon( logical(keep) );
    
    % save
    file  = [ datafolder_Eurogeog_country,'places_gecon.mat'];
    save(file,'places_gecon'); 
    
    %% plot
    %{
    h=figure;
    
    mapshow(places_gecon,'FaceColor','white')
    mapshow(country_bounds,'FaceColor','none')
    for i=1:length(places_gecon)
        text( places_gecon(i).Xcorner+1/2,places_gecon(i).Ycorner+1/2,...
            num2str( round(places_gecon(i).gdp*10)/10 ),...
            'color','red')
    end
    title(country) 
    if ismac
        saveas(h, [datafolder_gecon,'maps/',country,'_map_gecon.pdf']) 
    else
        saveas(h, [datafolder_gecon,'maps\',country,'_map_gecon.pdf']) 
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
    
    % load gecon
    file  = [ datafolder_Eurogeog_country,'places_gecon.mat'];
    temp = load(file); 
    
    % save temp file
    places_gecon_temp{country_n} = temp.places_gecon;
    for i=1:length(places_gecon_temp{country_n})
        places_gecon_temp{country_n}(i).country = country;
    end
    
end

i=1;
for country_n=1:Ncountries
    for j=1:length( places_gecon_temp{country_n} )
        places_gecon(i).Geometry =  places_gecon_temp{country_n}(j).Geometry;
        places_gecon(i).X =  places_gecon_temp{country_n}(j).X;
        places_gecon(i).Y =  places_gecon_temp{country_n}(j).Y;
        places_gecon(i).Xcorner =  places_gecon_temp{country_n}(j).Xcorner;
        places_gecon(i).Ycorner =  places_gecon_temp{country_n}(j).Ycorner;
        places_gecon(i).population =  places_gecon_temp{country_n}(j).population;
        places_gecon(i).gdp =  places_gecon_temp{country_n}(j).gdp;
        places_gecon(i).country_temp =  places_gecon_temp{country_n}(j).country;
        i=i+1;
    end
end

% combine cells that appear twice if they are on boundaries

% assign a unique id to each cell based on its geographic coordinates 
for i=1:length(places_gecon)
    XY(i,:) = [ places_gecon(i).Xcorner places_gecon(i).Ycorner ];
end      
[ uniqueXY,A,B ] = unique(XY,'rows');
for i=1:length(places_gecon)
    places_gecon(i).cell_id = B(i); 
end

% identify cells on frontier
for i=1:length(places_gecon)
    places_gecon(i).frontier = ( numel(find(B==B(i)))>1 );  % =1 if cell appears more than once, meaning it is on a frontier
end

% for each cell on frontier, assign the country in which this cell has the largest population
for i=1:length(places_gecon)
    if places_gecon(i).frontier == 1
        I = find( [places_gecon.cell_id]==places_gecon(i).cell_id );
        I = I([places_gecon(I).population]==max([places_gecon(I).population]));
        places_gecon(i).country = places_gecon(I).country_temp;
        if i==I
            places_gecon(i).keep=1;  % duplicated cell with highest population among duplicates
        else
            places_gecon(i).keep=0;
        end
    else   % non-duplicated cell
        places_gecon(i).country = places_gecon(i).country_temp;
        places_gecon(i).keep=1;
    end
end

% assign total income to each unique cell
for i=1:numel(A)
   
    inc_cell = sum( [ places_gecon( [places_gecon.cell_id]==i ).gdp ] );
    
    for j=find([places_gecon.cell_id]==i)
        places_gecon(j).tot_inc = inc_cell;
    end
    
end

% identify duplicated cells
duplicate_ind = setdiff( (1:length(places_gecon)),A );
places_gecon_duplicates = places_gecon(duplicate_ind);

% keep unique cells
places_gecon = places_gecon([places_gecon.keep]==1);

% clean up
for i=1:length(places_gecon)       
    places_gecon(i).gdp = places_gecon(i).tot_inc;
    places_gecon(i).rel_gdp = places_gecon(i).gdp/mean( [places_gecon.tot_inc ] );
end
places_gecon = rmfield( places_gecon,{'tot_inc','cell_id','keep','country_temp'});

% save
save([ datafolder_connected_countries,'places_gecon.mat'],'places_gecon');

% load connected country bounds
file = [ datafolder_connected_countries,'connected_countries_bound_many.mat'];
load(file);

% extra figures
%{
%duplicated cells
h=figure;   
mapshow(places_gecon,'FaceColor','red','EdgeColor','red')
mapshow(places_gecon([places_gecon.frontier]==1),'FaceColor','black','EdgeColor','black')
mapshow(country_bounds,'FaceColor','none','EdgeColor','white')

for i=1:length(places_gecon)
    if places_gecon(i).frontier==1
        text( places_gecon(i).Xcorner+1/2,places_gecon(i).Ycorner+1/2,places_gecon(i).country,...
              'Color','white' )
    end
end

% total income for unique cells
h=figure;
Colors = makesymbolspec('Polygon', {'gdp', ...
    [0 max( [places_gecon.gdp] )],...
    'FaceColor', jet(numel(places_gecon)) });
mapshow(places_gecon,'SymbolSpec',Colors) 
mapshow(country_bounds,'FaceColor','none','EdgeColor','white')    
saveas(h, ([ datafolder_connected_countries,'places_gecon.pdf'] ) )
%}