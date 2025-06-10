% Discretize map for European exercise

function [ code_control,places_grid,gridmap,edges,places_gecon ] = create_gridmap_connected( code_control,country_bounds,places_sedac,...
                                                                                             relief_etopo,places_gecon )
                                                                                                                                                              
% use 1-degree grid from G-econ as baseline grid
minX = min( [places_gecon.X] );
maxX = max( [places_gecon.X] );
minY = min( [places_gecon.Y] );
maxY = max( [places_gecon.Y] );

N_grid_X = ( maxX-minX );
N_grid_Y = ( maxY-minY );

% evenly spaced points over grid
Xgrid = linspace( minX,maxX,N_grid_X+1 );
Ygrid = linspace( minY,maxY,N_grid_Y+1 );
      
counter=1;
mat_counter = zeros( N_grid_X,N_grid_Y );
gridmap = [];
for i=1:N_grid_X
    for j=1:N_grid_Y
        gridmap(counter).Geometry = 'Polygon';
        gridmap(counter).X = [ Xgrid(i) Xgrid(i) Xgrid(i+1) Xgrid(i+1) NaN ];
        gridmap(counter).Y = [ Ygrid(j) Ygrid(j+1) Ygrid(j+1) Ygrid(j) NaN ];
        gridmap(counter).counter = counter;
        mat_counter(i,j) = counter;
        counter=counter+1;
    end    
end

%{
% show baseline grid and Gecon
mapshow(gridmap)
mapshow(places_gecon,'FaceColor','black')
%}
%%
%{
figure
mapshow(gridmap) 
mapshow( country_bounds,'FaceColor','Black' )
mapshow(relief_etopo)
%}

%% assign altitude and ruggedness to each cell

% X and Y coordinates in topo file
places_X_topo = cell2mat({relief_etopo.centerX}');
places_Y_topo = cell2mat({relief_etopo.centerY}');

% topo variables
av_alt = cell2mat( {relief_etopo.av_alt} );
sd_alt = cell2mat( {relief_etopo.sd_alt} );
rugged = cell2mat( {relief_etopo.rugged} );

for j=1:length( gridmap )

        % roughness
        in_j = inpolygon( places_X_topo, places_Y_topo, gridmap(j).X,gridmap(j).Y  ); % indicate whether each 0.25 degree cell in Etopo belongs in the cell
        
        gridmap(j).av_alt = mean( av_alt( in_j ) );  
        gridmap(j).sd_alt = mean( sd_alt( in_j ) );  
        gridmap(j).rugged = mean( rugged( in_j ) ); 

end
        
        
%% assign income and country to each cell from G-Econ

% assign quantiles
N_cat_income = 20;
quantiles_income = quantile( cell2mat({places_gecon.gdp}),N_cat_income );
for i=1:length(places_gecon)
    [ places_gecon(i).quantiles ] = sum( places_gecon(i).gdp>quantiles_income )+1;
end

% map cells in gridmap to cells in gecon
for i=1:length(gridmap)
    gridmap(i).Xcorner = gridmap(i).X(1);
    gridmap(i).Ycorner = gridmap(i).Y(1);
    gridmap(i).keep_gecon = 0;
end

% assign income and other variables from Gecon
for i=1:length(places_gecon)
    Ix = find( [gridmap.Xcorner]==places_gecon(i).Xcorner )';
    Iy = find( [gridmap.Ycorner]==places_gecon(i).Ycorner )';
    I = intersect(Ix,Iy);
    gridmap(I).income = places_gecon(i).gdp;
    gridmap(I).country = places_gecon(i).country;

    temp = name2code( places_gecon(i).country );
    gridmap(I).ccode = temp{1};
    gridmap(I).country_n = temp{2};

    gridmap(I).quantiles = places_gecon(i).quantiles;
    gridmap(I).frontier = places_gecon(i).frontier;
    gridmap(I).keep_gecon = 1;
end

for i=1:length(gridmap)
    if gridmap(i).keep_gecon == 0;
        gridmap(i).income = 0;
    end
end

gridmap = rmfield( gridmap,'keep_gecon' );

%{
figure
mapshow(gridmap( [gridmap.income]>0 ))
mapshow(places_gecon,'FaceColor','Black')
%}


%% determine neighbors
    counter=1;
    for i=1:N_grid_X

        for j=1:N_grid_Y

            if i==1 && j==1

                gridmap(counter).neighbors = [ mat_counter(i,j+1) mat_counter(i+1,j)  ];  % vertical-horizontal
                gridmap(counter).neighbors_D = mat_counter(i+1,j+1) ;                  % diagonal
                  
            elseif i==1 && j<N_grid_Y

                gridmap(counter).neighbors = [ mat_counter(i,j+1) mat_counter(i+1,j) mat_counter(i,j-1) ];
                gridmap(counter).neighbors_D = [ mat_counter(i+1,j+1) mat_counter(i+1,j-1) ];

            elseif i==1 && j==N_grid_Y

                gridmap(counter).neighbors = [ mat_counter(i+1,j) mat_counter(i,j-1)  ];
                gridmap(counter).neighbors_D = mat_counter(i+1,j-1);

            end

            if i>1 && i<N_grid_X && j==1

                gridmap(counter).neighbors = [ mat_counter(i-1,j) mat_counter(i,j+1) mat_counter(i+1,j)  ];
                gridmap(counter).neighbors_D = [ mat_counter(i-1,j+1) mat_counter(i+1,j+1)  ];

            elseif i>1 && i<N_grid_X && j<N_grid_Y

                gridmap(counter).neighbors = [ mat_counter(i-1,j) mat_counter(i,j+1) mat_counter(i+1,j) mat_counter(i,j-1) ];
                gridmap(counter).neighbors_D = [ mat_counter(i-1,j-1) mat_counter(i+1,j+1) mat_counter(i+1,j-1) mat_counter(i-1,j-1) ];

            elseif i>1 && i<N_grid_X && j==N_grid_Y

                gridmap(counter).neighbors = [ mat_counter(i+1,j) mat_counter(i,j-1) mat_counter(i-1,j) ];
                gridmap(counter).neighbors_D = [ mat_counter(i+1,j-1) mat_counter(i-1,j-1) ];

            end


            if i==N_grid_X && j==1

                gridmap(counter).neighbors = [ mat_counter(i-1,j) mat_counter(i,j+1)  ];
                gridmap(counter).neighbors_D = mat_counter(i-1,j+1 ) ;

            elseif i==N_grid_X && j<N_grid_Y

                gridmap(counter).neighbors = [ mat_counter(i,j-1) mat_counter(i-1,j) mat_counter(i,j+1)   ];
                gridmap(counter).neighbors_D = [ mat_counter(i-1,j-1) mat_counter(i-1,j+1) ];

            elseif i==N_grid_X && j==N_grid_Y

                gridmap(counter).neighbors = [ mat_counter(i,j-1) mat_counter(i-1,j) ];
                gridmap(counter).neighbors_D = mat_counter(i-1,j-1) ;

            end

            counter=counter+1;

        end

    end

%% centroids

% X and Y coordinates of each city
places_X = cell2mat({places_sedac.X}');
places_Y = cell2mat({places_sedac.Y}');

% population vector
pop_mat = cell2mat( {places_sedac.population}' ); % vector with city population

% identify cells in gridmap within the country map
keep = ones( N_grid_X*N_grid_Y,1 );
for j=1:N_grid_X*N_grid_Y
    temp_keep = ones( length(country_bounds),1 );
    for n=1:length(country_bounds)
        if sum( inpolygon( gridmap(j).X(1:4),gridmap(j).Y(1:4),country_bounds(n).X,country_bounds(n).Y ) ) == 0           
            temp_keep(n) = 0;  % if there is no intersection between square and bounds, drop the square        
        end
    %mapshow( gridmap(j),'FaceColor',(1-keep(j))*[ 1 0 0 ]+keep(j)*[ 1 1 1] ); % white=keep; red=drop
    %pause;
    end
    if max( temp_keep )==0
        keep(j) = 0;
    end
end

%{
close
mapshow(country_bounds)
mapshow(gridmap(logical(keep)),'FaceColor','none')
%}

%%
               
% recompute areas and population
for j=1:length( gridmap )
        
    if keep(j)==1
        
        % old coordinates of polygon
        x_old = gridmap(j).X;
        y_old = gridmap(j).Y;
        area_full_cell = areaint( y_old,x_old );
        
        % new coordinates of each polygon
        [x,y] = polybool( 'intersection',x_old,y_old,country_bounds.X,country_bounds.Y );    
        gridmap(j).X=x;
        gridmap(j).Y=y;
        
        % total population
        in_j = inpolygon( places_X, places_Y, gridmap(j).X,gridmap(j).Y  ); % indicator that each city is in cell j
        gridmap(j).population = in_j'*pop_mat;  
        
        % total area
        gridmap(j).area = sum( areaint( gridmap(j).Y,gridmap(j).X ) );
        
        % drop the cell to construct places structure?
        gridmap(j).drop = 0;            
        if ( gridmap(j).population==0 )|| ...                                % drop if no population 
           ( gridmap(j).area<=0.25*area_full_cell ) || ...            % drop if total area of cell <=25% of full cell
           ( max( areaint( gridmap(j).Y,gridmap(j).X ) )<=0.1*area_full_cell )     % drop if each connected polygon in cell <=10% of full cell
           keep(j) = 0;
        end

        % if cell remains, compute population centroid
        if keep(j)==1 % pick centroids of places with positive population
            
            if code_control.pick_largest   % pick largest city in each cell

                [ ~,j_largest ] = max( in_j.*pop_mat );  % index of largest city in cell j
                gridmap( j ).coordinates = [ places_X( j_largest ) places_Y( j_largest ) ];

            else  % pick centroids
                    
                gridmap( j ).coordinates = [ ( inpolygon( places_X, places_Y, gridmap(j).X,gridmap(j).Y  ).*places_X )'*pop_mat/gridmap(j).population ...
                                             ( inpolygon( places_X, places_Y, gridmap(j).X,gridmap(j).Y  ).*places_Y )'*pop_mat/gridmap(j).population  ];           

                % centroid can be outside the country for non-convex regions, if so pick largest city
                if ~inpolygon( gridmap(j).coordinates(1),gridmap(j).coordinates(2),country_bounds.X,country_bounds.Y )          
                        [ ~,j_largest ] = max( in_j.*pop_mat );  % index of largest city in cell j
                        gridmap( j ).coordinates = [ places_X( j_largest ) places_Y( j_largest ) ];   
                end

            end
            
        end
        
    end
    
end

% add index for whether node will be kept
for j=1:length(gridmap)
    gridmap(j).keep = keep(j);
end

% for each node, indicate whether some neighbor is kept
for j=1:length(gridmap)
    gridmap(j).keep2 = sum( keep( [ gridmap(j).neighbors gridmap(j).neighbors_D ] ) )>0;
end


%% adjust grid and edges to surviving neighbors

% keep surviving cells
surviving = find( keep );
gridmap = gridmap( surviving );   % grid with surviving nodes

% adjust set of neighbors
for j=1:length( gridmap )
    
    % new counter
    gridmap(j).newcounter = j;
    
    % exclude horizontal-vertical neighbors that were dropped when drawing intersection with boundaries
    temp = gridmap(j).neighbors;
    gridmap(j).neighbors = [];    
    for i=1:length( temp )
       
        if sum( temp(i)==surviving ) == 1  % neighbor is in surviving set
            
            gridmap(j).neighbors = [ gridmap(j).neighbors temp(i) ];
            
        end
        
    end    
    
    % exclude diagonal neighbors that were dropped when drawing intersection with boundaries
    if code_control.diagonals==1
        temp = gridmap(j).neighbors_D;
        gridmap(j).neighbors_D = [];    
        for i=1:length( temp )

            if sum( temp(i)==surviving ) == 1  % neighbor is in surviving set

                gridmap(j).neighbors_D = [ gridmap(j).neighbors_D temp(i) ];

            end


        end 
    end  
end

% relabel the neighbors using new_counter
temp = cell2mat({gridmap.counter})';
for j=1:length( gridmap )
    
    % horizontal-vertical
    for i=1:length( gridmap(j).neighbors )
        
        gridmap(j).neighbors(i) = find( temp==gridmap(j).neighbors(i),1,'first' );  % replace each neighbor for the current position of that neighbor in gridmap
        
    end
    
    % diagonal
    if code_control.diagonals==1
        for i=1:length( gridmap(j).neighbors_D )

            gridmap(j).neighbors_D(i) = find( temp==gridmap(j).neighbors_D(i),1,'first' );  % replace each neighbor for the current position of that neighbor in gridmap

        end

        gridmap(j).counter = gridmap(j).newcounter;
    end

end

% remove new counter as field
gridmap = rmfield(gridmap,'newcounter');

%% construct the edges

% vertical-horizontal grid

    % vectors with the origin and destination nodes
    gridmap_s = [];
    gridmap_t = [];
    for j=1:length(gridmap)
        gridmap_s = [ gridmap_s;gridmap(j).counter*ones( length( gridmap(j).neighbors ),1 ) ];
        gridmap_t = [ gridmap_t;gridmap(j).neighbors' ];
    end
    temp = [ gridmap_s gridmap_t];

    % keep unique horizontal-vertical edges
    for i=1:length( temp )    
        temp(i,:) = [ min( temp(i,:) ) max( temp(i,:) ) ];    
    end    
    [ unique_edges,~,~ ] = unique( temp,'rows');  % we will use unique_edges to build the vertical-horizontal graph

    %{
    %keep only edges that are fully contained in map
    keep = ones( length(unique_edges),1 );
    for j=1:length( unique_edges )
        
        init = gridmap( unique_edges(j,1) ).coordinates;    % starting coordinates in edge
        final = gridmap( unique_edges(j,2) ).coordinates;   % final coordinates in edge
        
        tempN = 10;
        tempX = linspace(init(1),final(1),tempN);
        tempY = linspace(init(2),final(2),tempN);
        
        keep(j) = ( sum( inpolygon( tempX,tempY,country_bounds.X,country_bounds.Y ) )==tempN );  % keep if entire edge is within  country
        
    end  
    unique_edges = unique_edges( logical(keep),: );
    %}
    
% diagonals    
if code_control.diagonals==1   

    % vectors with the origin and destination nodes
    gridmap_s = [];
    gridmap_t = [];
    for j=1:length(gridmap)
        gridmap_s = [ gridmap_s;gridmap(j).counter*ones( length( gridmap(j).neighbors_D ),1 ) ];
        gridmap_t = [ gridmap_t;gridmap(j).neighbors_D' ];
        temp = [ gridmap_s gridmap_t];
    end

    % keep unique diagonal edges
    for i=1:length( temp )    
        temp(i,:) = [ min( temp(i,:) ) max( temp(i,:) ) ];    
    end    
    [ unique_edges_D,~,~ ] = unique( temp,'rows');  % we will use unique_edges to build the diagonal graph
    
    %{
    % keep only edges fully contained in map
    keep = ones( length(unique_edges_D),1 );
    for j=1:length( unique_edges_D )
        
        init = gridmap( unique_edges_D(j,1) ).coordinates;    % starting coordinates in edge
        final = gridmap( unique_edges_D(j,2) ).coordinates;   % final coordinates in edge
        
        tempN = 10;
        tempX = linspace(init(1),final(1),tempN);
        tempY = linspace(init(2),final(2),tempN);
        
        keep(j) = ( sum( inpolygon( tempX,tempY,country_bounds.X,country_bounds.Y ) )==tempN );  % keep if entire edge is within country
        
    end  
    unique_edges_D = unique_edges_D( logical(keep),: );
    %}
end

% re-define all neighbors after eliminating edges
for j=1:length(gridmap)
    
    gridmap(j).neighbors = sort( [ unique_edges( unique_edges(:,1)==j,2 );...
                                   unique_edges( unique_edges(:,2)==j,1 );...
                                   unique_edges_D( unique_edges_D(:,1)==j,2 );...
                                   unique_edges_D( unique_edges_D(:,2)==j,1 ) ]' );
                               
    
end
gridmap = rmfield(gridmap,'neighbors_D');

%% idenfity cells on boundary
for i=1:length(gridmap)
   if length( gridmap(i).neighbors )<8
       gridmap(i).boundary = 1;
   else
       gridmap(i).boundary = 0;
   end
end

%% create the union of all components of gridmap
%{
X=country_bounds.X( ~isnan( [country_bounds.X] ) );
Y=country_bounds.Y( ~isnan( [country_bounds.Y] ) );

K=boundary(X,Y,1);

close
mapshow(X(K),Y(K))
mapshow(gridmap)
%}
%%
%{
close
figure
%mapshow( country_bounds,'FaceColor','Black' );
subplot(1,2,1)
mapshow( gridmap ,'FaceColor','white','EdgeColor','black' )
mapshow( gridmap([gridmap.boundary]==1),'FaceColor','red','EdgeColor','black' )
subplot(1,2,2)
mapshow( gridmap ,'FaceColor','white','EdgeColor','black' )
mapshow( gridmap([gridmap.frontier]==1),'FaceColor','red','EdgeColor','black' )
%}
%% show map and surviving cells
%{
figure
mapshow( country_bounds,'FaceColor','Black' );
mapshow(gridmap) 
%}
%% assign income per capita
for i=1:length(gridmap)
    
    % income per capita
    gridmap(i).y = gridmap(i).income/gridmap(i).population;  %income per capita
    if isnan(gridmap(i).y)
        gridmap(i).y = 0;
        gridmap(i).y = cast(gridmap(i).y,'single');
    end
                                       
end

% winsorize
ytemp = cell2mat( {gridmap.y} );
if length(gridmap)>20
    maxquantile = max( quantile(ytemp,20) );
elseif length(gridmap)>10
    maxquantile = max( quantile(ytemp,10) );
elseif length(gridmap)<=10
    maxquantile = max( quantile(ytemp,5) );
end

for i=1:length(gridmap)
    
    if gridmap(i).y>maxquantile
        gridmap(i).y = 0.99*maxquantile;
        gridmap(i).income = gridmap(i).y*gridmap(i).population;
    end 
    
end
%tabulate( cell2mat( {gridmap.y } ))


%% assign quantiles

% population
temp =  cell2mat({ gridmap.population })';
quantiles_pop_grid = quantile( temp( temp>0 ),code_control.N_cat );  % assign quantiles conditional on positive population in the square
for i=1:length(gridmap)
    if temp(i)>0
        [ gridmap(i).quantiles ] = sum( gridmap(i).population>quantiles_pop_grid )+1;
    elseif temp(i)==0
        [ gridmap(i).quantiles ] = 0;
    end
end

% income
temp =  cell2mat({ gridmap.y })';
quantiles_income_grid = quantile( temp( temp>0 ),code_control.N_cat );  % assign quantiles conditional on positive population in the square
for i=1:length(gridmap)
    if temp(i)>0
        [ gridmap(i).quantiles_inc ] = sum( gridmap(i).y>quantiles_income_grid )+1;
    elseif temp(i)==0
        [ gridmap(i).quantiles_inc ] = 0;
    end
end   

% assign altitude and ruggedness
temp1 =  cell2mat({ gridmap.av_alt })';
temp2 =  cell2mat({ gridmap.rugged })';
for i=1:length(gridmap)
    [ gridmap(i).quantiles_av_alt ] = sum( gridmap(i).av_alt>quantile( temp1,code_control.N_cat ) )+1;
    [ gridmap(i).quantiles_rugged ] = sum( gridmap(i).rugged>quantile( temp2,code_control.N_cat ) )+1;
end

%tabulate(cell2mat({gridmap.quantiles})); % check out frequency of population quantiles

%% create a new places file from gridmap
for j=1:length(gridmap);
    
    places_grid(j).Geometry = 'Point';
    places_grid(j).X = gridmap(j).coordinates(1);
    places_grid(j).Y = gridmap(j).coordinates(2);
    places_grid(j).population = gridmap(j).population;
    places_grid(j).income = gridmap(j).income;
    places_grid(j).income = cast(places_grid(j).income,'single');
    places_grid(j).y = gridmap(j).y;
    %places_grid(j).quantiles = gridmap(j).quantiles;  %pop_quantile
    places_grid(j).neighbors = gridmap(j).neighbors;
    places_grid(j).av_alt = gridmap(j).av_alt;
    places_grid(j).sd_alt = gridmap(j).sd_alt;
    places_grid(j).rugged = gridmap(j).rugged;
    places_grid(j).country = gridmap(j).country;
    places_grid(j).country_n = gridmap(j).country_n;
    places_grid(j).ccode = gridmap(j).ccode;
    
end

edges.unique_edges = unique_edges;
edges.unique_edges_D = unique_edges_D;

%% plot income
%{
close
subplot(2,2,1)
gdp = cell2mat({places_gecon.gdp});
mapshow( places_gecon,'FaceColor','White')
for j=1:length(places_gecon)
    rel = gdp(j)/max(gdp);
    mapshow( places_gecon(j),'FaceColor',[0 0 1]*rel)
end
title('G-Econ (Raw Income Data)')
                margin = 0.5;
                Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
                Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
                set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )


subplot(2,2,2)
gdp = cell2mat({places_grid.income});
mapshow( gridmap,'FaceColor','White')
for j=1:length(places_grid)
    rel = gdp(j)/max(gdp);
    mapshow( gridmap(j),'FaceColor',[0 0 1]*rel)
end
                margin = 0.5;
                Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
                Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
                set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )
title('Income in our Data')               
                
                
subplot(2,2,3)
pop = cell2mat({places_grid.population});
mapshow( gridmap,'FaceColor','White')
for j=1:length(places_grid)
    rel = pop(j)/max(pop);
    mapshow( gridmap(j),'FaceColor',[0 0 1]*rel)
end
                margin = 0.5;
                Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
                Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
                set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )
title('Population in our Data')

subplot(2,2,4)
gdppc = cell2mat({places_grid.y});
mapshow( gridmap,'FaceColor','White')
for j=1:length(places_grid)
    rel = gdppc(j)/max(gdppc);
    mapshow( gridmap(j),'FaceColor',[0 0 1]*rel)
end
[~,i] = max(gdppc)
pop(i)
mapshow( gridmap(i),'FaceColor',[1 0 0] )

                margin = 0.5;
                Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
                Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
                set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )
title('Income per capita in our Data')

%}

