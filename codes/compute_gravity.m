function gravity_coeff = compute_gravity( bilateral_trade_mat,bilateral_dist ) 

if size(bilateral_trade_mat,1)~=size(bilateral_trade_mat,2) || size(bilateral_dist,1)~=size(bilateral_dist,2)
    error('trade matrix not square')
end

if size(bilateral_trade_mat)~=size(bilateral_dist)
    error('Trade and Distance Matrixes have different sizes')
end

% number of NUTS
nb_NUTS = size(bilateral_dist,1);

% long
bilateral_trade_long = zeros( nb_NUTS*nb_NUTS,5+2*nb_NUTS );  % o,d,log(flow),log(distance),...
                                                              % constant
                                                              % {FEs_origin},{FEs_destination}
                                                             
keep = ones( nb_NUTS*nb_NUTS,1 );    %which rows to keep  

j=1;
for o=1:nb_NUTS
    for d=1:nb_NUTS
        bilateral_trade_long( j,1:5 ) = [ o d log( bilateral_trade_mat(o,d) ) log( bilateral_dist(o,d) ) 1 ];
        bilateral_trade_long( j,4+o ) = 1;  % FE origin
        bilateral_trade_long( j,4+nb_NUTS+d ) = 1;   % FE destination
        if (o==d) || bilateral_trade_mat(o,d)==0
            keep(j) = 0;
        end
        j=j+1;
    end
end

% keep non-zero trade flows and exclude trade with own region
bilateral_trade_long = bilateral_trade_long( logical(keep),: );

% run gravity with fixed effects
X = bilateral_trade_long( :,4:size(bilateral_trade_long,2) );
y = bilateral_trade_long( :,3 );

% save the coefficients
mdl=fitlm(X,y)
gravity_coeff.b = mdl.Coefficients.Estimate(2);
gravity_coeff.se = mdl.Coefficients.SE(2);