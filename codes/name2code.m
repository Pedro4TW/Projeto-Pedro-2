function code = name2code( country )

    if strcmp( country,'Austria' );
        code = {'AT',1};
    elseif strcmp( country,'Belgium' )
        code = {'BE',2};           
    elseif strcmp( country,'Czech Republic' );
        code = {'CZ',3};
    elseif strcmp( country,'Germany' );
        code = {'DE',4};
    elseif strcmp( country,'Denmark' );
        code = {'DK',5};
    elseif strcmp( country,'Spain' );
        code = {'ES',6};
    elseif strcmp( country,'France' );
        code = {'FR',7};  
    elseif strcmp( country,'Hungary' );
        code = {'HU',8};
    elseif strcmp( country,'Italy' );
        code = {'IT',9};
    elseif strcmp( country,'Netherlands' );
        code = {'NL',10};
    elseif strcmp( country,'Portugal' );
        code = {'PT',11};
    elseif strcmp( country,'Slovenia' );
        code = {'SI',12};
    elseif strcmp( country,'Slovakia' );
        code = {'SK',13};
    elseif strcmp( country,'Switzerland'  );
        code = {'CHLI',14};
    end