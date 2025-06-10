temp = dir(datafolder_Eurogeog);

for k = length(temp):-1:1 % this method should work for Mac/Linux/Windows
    if ~temp(k).isdir || temp(k).name(1) == '.' % remove non-folders and hidden folders
        temp(k) = [];
    end
end

countries = {temp.name}';

% keep large countries
keep_countries = {'ES','FR'};

keep = zeros(length(keep_countries),1);
for j=1:length(keep)
    keep(j) = find( strcmpi( countries,keep_countries(j) ) );
end

countries = countries(keep);

% add names
Ncountries = length(countries);
country_names = cell( Ncountries,3 );
for j=1:Ncountries
    country_names{j,1} = j;
    country_names{j,2} = char( countries(j) );
    country_names{j,3}  = icc2name( char( countries(j) ) );
end

country_names 