temp = dir(datafolder_Eurogeog);

for k = length(temp):-1:1 % this method should work for Mac/Linux/Windows
    if ~temp(k).isdir || temp(k).name(1) == '.' % remove non-folders and hidden folders
        temp(k) = [];
    end
end

countries_connected = {temp.name}';

% exclude from the analysis the countries with missing lane data in EGM
exclude = {'KS','MT','BG','GB','GR','HR','NO','RO','UA','IS','PL','SE','EE','RS',...
           'CY','FI','GE','IE','LT','LV','MK','ND','MD'};
       
       % in the first run, we excluded: 'SI','HU','SK','CZ'

for j=1:length(exclude)
    countries_connected = countries_connected( ~strcmp( countries_connected,exclude(j) ) );
end

% add names
Ncountries = length(countries_connected);
country_names_connected = cell( Ncountries,3 );
for j=1:Ncountries
    country_names_connected{j,1} = j;
    country_names_connected{j,2} = char( countries_connected(j) );
    country_names_connected{j,3}  = icc2name( char( countries_connected(j) ) );
end

country_names_connected 
