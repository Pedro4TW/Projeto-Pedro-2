% keep rods within country bounds

function roads = bound_roads( roads,country_bounds,use_hull )

if use_hull
    
    X=country_bounds.X( ~isnan( [country_bounds.X] ) );
    Y=country_bounds.Y( ~isnan( [country_bounds.Y] ) );        

    K=boundary(X,Y,0.8);
    X = X(K)';
    Y = Y(K)';    
    
else

    X = country_bounds.X';
    Y = country_bounds.Y';

end
    
for j=1:length( roads )    
    
    [ j length(roads) ]
    
    tempX = roads(j).X( 1:( length( roads(j).X )-1 ) );
    tempY = roads(j).Y( 1:( length( roads(j).Y )-1 ) );        
    
    keep = inpolygon( tempX,tempY,X,Y );
    
    roads(j).X = [ tempX( keep ) NaN];
    roads(j).Y = [ tempY( keep ) NaN];

end

