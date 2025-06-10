temp = dir(datafolder_Eurogeog);

for k = length(temp):-1:1 % this method should work for Mac/Linux/Windows
    if ~temp(k).isdir || temp(k).name(1) == '.' % remove non-folders and hidden folders
        temp(k) = [];
    end
end

countries = {temp.name}';

% exclude from the analysis the countries with missing lane data and some
% additional countries
exclude = {'KS','MT','BG','GB','GR','HR','NO','RO','UA','IS','PL','SE','EE','RS',...
           'GE','MD','ND','CY','LU',... % countries without NUTS classification
           'LT','LV','MK','SI'};        % further exclude countries with just 1 or 2 NUTS2

for j=1:length(exclude)
    countries = countries( ~strcmp( countries,exclude(j) ) );
end

% add names
Ncountries = length(countries);
country_names = cell( Ncountries,3 );
for j=1:Ncountries
    country_names{j,1} = j;
    country_names{j,2} = char( countries(j) );
    country_names{j,3}  = icc2name( char( countries(j) ) );
end

country_names 