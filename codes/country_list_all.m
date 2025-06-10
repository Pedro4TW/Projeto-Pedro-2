temp = dir(datafolder_Eurogeog);

for k = length(temp):-1:1 % this method should work for Mac/Linux/Windows
    if ~temp(k).isdir || temp(k).name(1) == '.' % remove non-folders and hidden folders
        temp(k) = [];
    end
end

countries_all = {temp.name}';

% exclude from the analysis the countries with missing lane data in EGM
%exclude = {'KS','MT','BG','GB','GR','HR','NO','RO','UA','IS','PL','SE','EE','RS'};

% add names
Ncountries = length(countries_all);
country_names_all = cell( Ncountries,3 );
for j=1:Ncountries
    country_names_all{j,1} = j;
    country_names_all{j,2} = char( countries_all(j) );
    country_names_all{j,3}  = icc2name( char( countries_all(j) ) );
end

country_names_all 