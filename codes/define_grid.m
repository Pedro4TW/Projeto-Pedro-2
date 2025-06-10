% construct evenly spaced grid from boundary of a map

function [ minX,maxX,minY,maxY,N_grid_X,N_grid_Y ] = define_grid( country_bounds,cellsize )

% chose grid corners to the nearest integer

for i=1:length(country_bounds)
    minX(i) = ( country_bounds(i).BoundingBox( 1,1 )-floor( country_bounds(i).BoundingBox( 1,1 ) )>cellsize )*( floor( country_bounds(i).BoundingBox( 1,1 ) )+cellsize )+...
           ( country_bounds(i).BoundingBox( 1,1 )-floor( country_bounds(i).BoundingBox( 1,1 ) )<=cellsize )*floor( country_bounds(i).BoundingBox( 1,1 ) );
    minY(i) = ( country_bounds(i).BoundingBox( 1,2 )-floor( country_bounds(i).BoundingBox( 1,2 ) )>cellsize )*( floor( country_bounds(i).BoundingBox( 1,2 ) )+cellsize )+...
           ( country_bounds(i).BoundingBox( 1,2 )-floor( country_bounds(i).BoundingBox( 1,2 ) )<=cellsize )*floor( country_bounds(i).BoundingBox( 1,2 ) );
    maxX(i) = ( ceil( country_bounds(i).BoundingBox( 2,1 ) )-country_bounds(i).BoundingBox( 2,1 )>cellsize )*(  ceil( country_bounds(i).BoundingBox( 2,1 ) )-cellsize )+...
           ( ceil( country_bounds(i).BoundingBox( 2,1 ) )-country_bounds(i).BoundingBox( 2,1 )<=cellsize )*ceil( country_bounds(i).BoundingBox( 2,1 ) );
    maxY(i) = ( ceil( country_bounds(i).BoundingBox( 2,2 ) )-country_bounds(i).BoundingBox( 2,2 )>cellsize )*(  ceil( country_bounds(i).BoundingBox( 2,2 ) )-cellsize )+...
           ( ceil( country_bounds(i).BoundingBox( 2,2 ) )-country_bounds(i).BoundingBox( 2,2 )<=cellsize )*ceil( country_bounds(i).BoundingBox( 2,2 ) );
end

minX = min(minX);
minY = min(minY);
maxX = max(maxX);
maxY = max(maxY);

% grid
N_grid_X = ( maxX-minX )/cellsize;
N_grid_Y = ( maxY-minY )/cellsize;
