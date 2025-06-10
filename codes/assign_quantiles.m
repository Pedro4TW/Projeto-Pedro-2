function out = assign_quantiles( in )

temp = quantile( in,9 );  % assign quantiles conditional on positive population in the square
out = zeros( length( in ),1 );
for i=1:length( out )    
    out( i ) = sum( in(i)>temp )+1;    
end