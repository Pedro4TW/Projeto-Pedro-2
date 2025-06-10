function taxes_struct = compute_taxes (param,graph,res)

% -------------
% Compute taxes
% -------------

if param.cong==false % good-specific congestion
    Pjn_mat = permute(repmat(res.Pjn,[1 1 graph.J]),[1 3 2]);
    Pkn_mat = permute(repmat(res.Pjn,[1 1 graph.J]),[3 1 2]);
    
    num = repmat( param.beta.*graph.delta_tau.*res.Ijk.^(-param.gamma), [1 1 param.N])...
        .* res.Qjkn.^param.beta; 
    
    denom = 1 + repmat( (1+param.beta).*graph.delta_tau.*res.Ijk.^(-param.gamma), [1 1 param.N])...
        .* res.Qjkn.^param.beta; 
        
    tjkn= num ./ denom;
    PtjknQ = Pkn_mat .* tjkn .* res.Qjkn;
        
    % Make sure to set non-adjacent places to 0
    PtjknQ( repmat(~graph.adjacency,[1 1 param.N]) ) = 0;
    tjkn( repmat(~graph.adjacency,[1 1 param.N]) ) = 0;
else % cross-good congestion
    Pjn_mat = permute(repmat(res.Pjn,[1 1 graph.J]),[1 3 2]);
    Pkn_mat = permute(repmat(res.Pjn,[1 1 graph.J]),[3 1 2]);
    PCj_mat = repmat(res.PCj,[1 graph.J param.N]);
    
    Qjk = sum( permute(repmat( param.m, [1 graph.J graph.J]),[2 3 1]) .* res.Qjkn, 3);
    
    num = repmat( param.beta.*graph.delta_tau.*res.Ijk.^(-param.gamma), [1 1 param.N])...
        .* permute(repmat( param.m, [1 graph.J graph.J]),[2 3 1])...
        .* PCj_mat... 
        .* repmat( Qjk, [1 1 param.N] ).^param.beta; 
    
    denom = Pjn_mat + repmat( (1+param.beta).*graph.delta_tau.*res.Ijk.^(-param.gamma), [1 1 param.N])...
        .* permute(repmat( param.m, [1 graph.J graph.J]),[2 3 1])...
        .* PCj_mat... 
        .* repmat( Qjk, [1 1 param.N] ).^param.beta; 
    
    tjkn= num ./ denom; % compute ad-valorem taxes per unit sent in terms of good at destination
    PtjknQ = Pkn_mat .* tjkn .* res.Qjkn; % compute total tax revenue
        
    % Make sure to set non-adjacent places to 0
    PtjknQ( repmat(~graph.adjacency,[1 1 param.N]) ) = 0;
    tjkn( repmat(~graph.adjacency,[1 1 param.N]) ) = 0;
    
end
    
% -----------------
% Report statistics
% -----------------

GDPj = res.PCj.*res.Cj/param.alpha; % total GDP (including transfers) = PC*C+PH*H
value_added_tradeablej = sum(res.Pjn.*res.Yjn,2);

taxes_struct.tjkn = tjkn;
taxes_struct.PtjknQ = PtjknQ;
taxes_struct.mean_ad_valorem = mean ( tjkn(repmat(graph.adjacency==1,[1 1 param.N])) );
taxes_struct.std_ad_valorem = std ( tjkn(repmat(graph.adjacency==1,[1 1 param.N])) );
taxes_struct.max_ad_valorem = max ( tjkn(repmat(graph.adjacency==1,[1 1 param.N])) );
taxes_struct.min_ad_valorem = min ( tjkn(repmat(graph.adjacency==1,[1 1 param.N])) );

taxes_struct.total_t_over_GDP = mean( sum(sum(PtjknQ, 3),2) ./ GDPj );
taxes_struct.total_t_over_tradeable = mean( sum(sum(PtjknQ, 3),2) ./ value_added_tradeablej );

% fprintf('--------------------------------\n');
% fprintf('Statistics for ad-valorem taxes:\n');
% fprintf('Mean=%.3f, std=%.3f, min=%.3f, max=%.3f\n',mean_t,std_t,min_t,max_t);
% fprintf('Tax revenue as share of GDP=%.3f%%\n',mean_t_GDP);
% fprintf('Tax revenue as share of value added in tradeable=%.3f%%\n',mean_t_tradeable);

