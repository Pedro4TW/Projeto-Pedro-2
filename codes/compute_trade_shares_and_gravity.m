function [ intra_reg_trade_share_differentiated,gravity_coeff,bilateral_trade_mat ] = compute_trade_shares_and_gravity(param,g,results,homogeneous_index,bilateral_dist)

%% --------------------
% DIFFERENTIATED GOODS
% -------------------- 

% set of differentiated goods
diff_goods = setdiff( (1:param.N)',homogeneous_index ); % get the set of differentiated goods

%[NUTS_list,~,g.region_id] = unique(g.region); % define a new region index between 1 and nb_NUTS so that NUTS_list(g.region_id)=g.region
nb_NUTS = size(bilateral_dist,1); %length(NUTS_list);

% Create a vector that identifies which good is produced in each NUTS
good_to_nuts=zeros(param.N-1,1);
for n=1:param.N-1
    node_n = param.Zjn(:,diff_goods(n))~=0;  % node that produces the differentiated good diff_goods(n)
    good_to_nuts(n)=g.region( node_n );      % region that produces the differentiated good %g.region_id( node_n );  
end

% Recover trade data
if param.cong==true
    Cjn=results.Djn;
else
    Cjn=results.Cjn;
end

% Construct trade matrix
bilateral_trade_mat=zeros( nb_NUTS,nb_NUTS );
%{
for j=1:g.J
    for n=1:param.N-1
        bilateral_trade_mat( good_to_nuts(n),g.region_id(j) ) = bilateral_trade_mat( good_to_nuts(n),g.region_id(j) )...
                                                                +results.Pjn(j,diff_goods(n) )*Cjn(j,diff_goods(n) );
    end    
end
%}
for n=1:param.N-1
    
    for j=1:g.J
    
        node_n = param.Zjn(:,diff_goods(n))~=0;  % node that produces the differentiated good diff_goods(n)
        
        bilateral_trade_mat( good_to_nuts(n),g.region(j) ) = bilateral_trade_mat( good_to_nuts(n),g.region(j) )...
                                                                +results.Pjn(j,diff_goods(n) )*Cjn(j,diff_goods(n) );
                                                            
        bilateral_trade_mat_orig_price( good_to_nuts(n),g.region(j) ) = bilateral_trade_mat( good_to_nuts(n),g.region(j) )...
                                                                +results.Pjn(node_n,diff_goods(n) )*Cjn(j,diff_goods(n) );                                                    
    end   
    
end

% compute gravity coefficient
gravity_coeff = compute_gravity( bilateral_trade_mat,bilateral_dist )

% compute internal share
intra_reg_trade_share_differentiated=sum(diag(bilateral_trade_mat))/sum(bilateral_trade_mat(:))


intra_reg_trade_share_differentiated_orig_price=sum(diag(bilateral_trade_mat_orig_price))/sum(bilateral_trade_mat_orig_price(:))

%% bring homogeneous good to the calculation

% homogeneous consumption
P0C0 = results.Pjn(:,1).*Cjn(:,1);

% homogeneous production
P0Y0 = results.Pjn(:,1).*results.Yjn(:,1);

% intra-region flows
for nuts=1:length(unique(g.region))
    
       intra0(nuts) = min( sum( P0C0(g.region==nuts) ),sum( P0Y0(g.region==nuts) ) );
end

% intra-region trade with hom good
intra_reg_trade_share_diff_hom = ...
       ( sum( diag(bilateral_trade_mat_orig_price) )+sum( intra0 ) )/... 
       ( sum( bilateral_trade_mat_orig_price(:) )+sum( P0Y0 ) );

% intra-region trade hom good only
intra_reg_trade_share_hom = sum( intra0 )/sum( P0Y0 );
