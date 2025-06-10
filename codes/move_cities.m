% move the cenotroids to the road network

function places2 = move_cities( gridmap,places_grid,unique_nodes,unique_nodes_nat )
 
%%
places2 = places_grid;
    
    for j=1:length( places2 )
                
        % determine the set of national roads within the cell
        
        if isempty( unique_nodes_nat )
            
            nat_roads_in_cell = 0;
        
        else
            
            nat_roads_in_cell = inpolygon( unique_nodes_nat(:,1),unique_nodes_nat(:,2),...
                                       gridmap(j).X,gridmap(j).Y );  
                                   
        end
                                   
                              
        if sum( nat_roads_in_cell )>0    
            
            % closest point in the road network within cell of city 
            coords_nat_roads_in_cell = unique_nodes_nat( nat_roads_in_cell,: ); 
            new_coords_nodes = dsearchn( coords_nat_roads_in_cell,...
                                         [ places_grid(j).X places_grid(j).Y ] );  

            % move the city
            new_coords = coords_nat_roads_in_cell( new_coords_nodes,: );

            places2(j).X = new_coords(1);        % places2 is same as places_grid, but the centroid coordinates are on the transport network
            places2(j).Y = new_coords(2);
            
            % node on the road network corresponding to this place
            nodes_roads_in_cell = find( unique_nodes( :,1 ) == new_coords(1) );
            nodes_roads_in_cell = nodes_roads_in_cell( unique_nodes( nodes_roads_in_cell,2 ) == new_coords(2) );            
            places2(j).node = nodes_roads_in_cell;  
            
            places2(j).in_road_network = 1;      % =1 indicates that the place is on the road network                       

        else   % try all roads
            
            % determine the set of all roads within the cell
            roads_in_cell = inpolygon( unique_nodes(:,1),unique_nodes(:,2),...
                                       gridmap(j).X,gridmap(j).Y );

            coords_roads_in_cell =  unique_nodes( roads_in_cell,: );    
            nodes_roads_in_cell =  find( roads_in_cell==1 );
            
            if sum( roads_in_cell )>0    
            
                % closest point in the road network within cell of city 
                new_coords_nodes = dsearchn( coords_roads_in_cell,...
                                             [ places_grid(j).X places_grid(j).Y ] );  

                % move the city
                new_coords = coords_roads_in_cell( new_coords_nodes,: );

                places2(j).X = new_coords(1);                               % places2 is exactly like places_grid, but the centroid coordinates are on the transport network
                places2(j).Y = new_coords(2);
                places2(j).node = nodes_roads_in_cell( new_coords_nodes );  % node on the road network corresponding to this place
                places2(j).in_road_network = 1;                             % =1 indicates that the place is on the road network
            
            else
            
                places2(j).in_road_network = 0;  % place is not on the road network
                
            end
            
        end

    end 
    
    % convert to double
    for j=1:length(places2)
        places2(j).X = double( places2(j).X );
        places2(j).Y = double( places2(j).Y );
    end