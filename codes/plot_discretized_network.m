                          % boundaries
                          mapshow( country_bounds,...
                          'FaceColor','Black',...
                          'EdgeColor','White','LineWidth', 0.5 )                                                                                                                                            
  

                          % edges
                          for j=1:length( edges )

                                 if edges(j)~=0

                                        colorweight = abs( edges(j) )/max_bright;

                                        if edges(j)>0
                                            color =  min( bright_line+( 1-bright_line )*colorweight*color_growth,1 );  
                                        elseif edges(j)<0
                                            color =  min( bright_line+( 1-bright_line )*colorweight*color_shrink,1 ); 
                                        end

                                    mapshow( discretized_roads( j ),...
                                            'Color',min( color,1 ),...
                                            'LineWidth',abs( edges(j) )*adj_width )

                                    hold on;

                                end

                          end 
                             
                          % nodes
                          x = nodes;
                          for i=1:length( gridmap )
                              
                              if x(i)~=0

                                    norm_neg = 1/abs( mean(x(x<0)) )*1/2.5;   %gL(gL<0)*norm_neg 
                                    norm_pos = 1/max( x(x>0) )*2.5;           %gL(gL>0)*norm_neg 

                                    relsize = x(i);
                                    bring_fire = 0;

                                    color = min( max( -relsize*norm_neg,0 )*color_shrink+...
                                                      ( max( -relsize*norm_neg,0.5 )-0.5 )*bring_fire*[ 0 1 0 ],1 )+...
                                                 min( max( relsize*norm_pos,0 )*color_growth+...
                                                      ( max( relsize*norm_pos,0.5 )-0.5 )*bring_fire*[ 0 1 0 ],1 );
            
                                    mapshow( places2( i ),...
                                        'Marker','o',...
                                        'MarkerFaceColor',color,...
                                        'MarkerEdgeColor',[ 1 1 1 ],...
                                        'MarkerSize',markersize )  
                                    
                              end

                            end                               