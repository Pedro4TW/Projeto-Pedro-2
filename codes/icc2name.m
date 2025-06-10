function country = icc2name( country_icc )

    % country name for natural earth
    if strcmp( country_icc,'AT' )
        country = 'Austria';
    elseif strcmp( country_icc,'BE' )
        country = 'Belgium';            % Belgium not in NE
    elseif strcmp( country_icc,'BG' )
        country = 'Bulgaria';
    elseif strcmp( country_icc,'CHLI' )
        country = 'Switzerland';
    elseif strcmp( country_icc,'CY' )
        country = 'Cyprus';
    elseif strcmp( country_icc,'CZ' )
        country = 'Czech Republic';
    elseif strcmp( country_icc,'DE' )
        country = 'Germany';
    elseif strcmp( country_icc,'DK' )
        country = 'Denmark';
    elseif strcmp( country_icc,'EE' )
        country = 'Estonia';
    elseif strcmp( country_icc,'ES' )
        country = 'Spain';
    elseif strcmp( country_icc,'FI' )
        country = 'Finland';
    elseif strcmp( country_icc,'FR' )
        country = 'France';
    elseif strcmp( country_icc,'GB' )
        country = 'Great Britain';    
    elseif strcmp( country_icc,'GE' )
        country = 'Georgia';
    elseif strcmp( country_icc,'GR' )
        country = 'Greece';
    elseif strcmp( country_icc,'HR' )
        country = 'Croatia';
    elseif strcmp( country_icc,'HU' )
        country = 'Hungary';
    elseif strcmp( country_icc,'IE' )
        country = 'Ireland';
    elseif strcmp( country_icc,'IS' )
        country = 'Iceland';
    elseif strcmp( country_icc,'IT' )
        country = 'Italy';
    elseif strcmp( country_icc,'KS' )
        country = 'Kosovo';                  
    elseif strcmp( country_icc,'LT' )
        country = 'Lithuania';
    elseif strcmp( country_icc,'LU' )
        country = 'Luxembourg';
    elseif strcmp( country_icc,'LV' )
        country = 'Latvia';
    elseif strcmp( country_icc,'MD' )
        country = 'Moldova';
    elseif strcmp( country_icc,'MK' )  
        country = 'Macedonia';
    elseif strcmp( country_icc,'MT' )
        country = 'Malta';
    elseif strcmp( country_icc,'ND' )
        country = 'Northern Ireland';   
    elseif strcmp( country_icc,'NL' )
        country = 'Netherlands';
    elseif strcmp( country_icc,'NO' )
        country = 'Norway';
    elseif strcmp( country_icc,'PL' )
        country = 'Poland';
    elseif strcmp( country_icc,'PT' )
        country = 'Portugal';
    elseif strcmp( country_icc,'RO' )
        country = 'Romania';
    elseif strcmp( country_icc,'RS' )
        country = 'Serbia';
    elseif strcmp( country_icc,'SE' )
        country = 'Sweden';
    elseif strcmp( country_icc,'SI' )
        country = 'Slovenia';
    elseif strcmp( country_icc,'SK' )
        country = 'Slovakia';
    elseif strcmp( country_icc,'UA' )
        country = 'Ukraine';
    elseif strcmp( country_icc,'ALL' )
        country = 'conn_countries';
    end