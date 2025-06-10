function country = icc2egm( country_icc )

    % ICC in EGM
    if strcmp( country_icc,'AT' )
        country = 'AT';
    elseif strcmp( country_icc,'BE' )
        country = 'BE';            % Belgium not in NE
    elseif strcmp( country_icc,'BG' )
        country = 'BG';
    elseif strcmp( country_icc,'CHLI' )
        country = 'CH';
    elseif strcmp( country_icc,'CY' )
        country = 'CY';
    elseif strcmp( country_icc,'CZ' )
        country = 'CZ';
    elseif strcmp( country_icc,'DE' )
        country = 'DE';
    elseif strcmp( country_icc,'DK' )
        country = 'DK';
    elseif strcmp( country_icc,'EE' )
        country = 'EE';
    elseif strcmp( country_icc,'ES' )
        country = 'ES';
    elseif strcmp( country_icc,'FI' )
        country = 'FI';
    elseif strcmp( country_icc,'FR' )
        country = 'FR';
    elseif strcmp( country_icc,'GB' )
        country = 'GB';    
    elseif strcmp( country_icc,'GE' )
        country = 'GE';
    elseif strcmp( country_icc,'GR' )
        country = 'GR';
    elseif strcmp( country_icc,'HR' )
        country = 'HR';
    elseif strcmp( country_icc,'HU' )
        country = 'HU';
    elseif strcmp( country_icc,'IE' )
        country = 'IE';
    elseif strcmp( country_icc,'IS' )
        country = 'IS';
    elseif strcmp( country_icc,'IT' )
        country = 'IT';
    elseif strcmp( country_icc,'KS' )
        country = 'XK';                  
    elseif strcmp( country_icc,'LT' )
        country = 'LT';
    elseif strcmp( country_icc,'LU' )
        country = 'LU';
    elseif strcmp( country_icc,'LV' )
        country = 'LV';
    elseif strcmp( country_icc,'MD' )
        country = 'MD';
    elseif strcmp( country_icc,'MK' )  
        country = 'MK';
    elseif strcmp( country_icc,'MT' )
        country = 'MT';
    elseif strcmp( country_icc,'ND' )
        country = 'ND';   
    elseif strcmp( country_icc,'NL' )
        country = 'NL';
    elseif strcmp( country_icc,'NO' )
        country = 'NO';
    elseif strcmp( country_icc,'PL' )
        country = 'PL';
    elseif strcmp( country_icc,'PT' )
        country = 'PT';
    elseif strcmp( country_icc,'RO' )
        country = 'RO';
    elseif strcmp( country_icc,'RS' )
        country = 'RS';
    elseif strcmp( country_icc,'SE' )
        country = 'SE';
    elseif strcmp( country_icc,'SI' )
        country = 'SI';
    elseif strcmp( country_icc,'SK' )
        country = 'SK';
    elseif strcmp( country_icc,'UA' )
        country = 'UA';
    end