function [X,M] = recover_X_M( results,g )

Qjkn = results.Qjkn;
ngoods = size(Qjkn,3);
temp = zeros( g.J );
for j=1:ngoods
    temp = temp+Qjkn( :,:,j ).*repmat( results.Pjn( :,j ),1,g.J );
end
X = sum( temp,2 );  % total exports by origin

temp = zeros( g.J );
for j=1:ngoods
    temp = temp+Qjkn( :,:,j ).*repmat( results.Pjn( :,j )',g.J,1 );
end
M = sum( temp,1 )';  % total imports by origin