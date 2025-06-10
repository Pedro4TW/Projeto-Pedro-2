function [intra_reg_trade_share,intra_reg_trade_share_differentiated] = compute_trade_shares(param,g,results,homogeneous_index)

% --------------------
% DIFFERENTIATED GOODS
% --------------------

diff_goods = setdiff((1:param.N)',homogeneous_index); % get the set of differentiated goods

[NUTS_list,~,g.region_id] = unique(g.region); % define a new region index between 1 and nb_NUTS so that NUTS_list(g.region_id)=g.region
nb_NUTS = length(NUTS_list);

% Create a vector that identifies which good is produced in each NUTS
good_to_nuts=zeros(param.N-1,1);
for n=1:param.N-1
    list=find(param.Zjn(:,diff_goods(n))~=0);
    good_to_nuts(n)=g.region_id(list(1));
end

% Recover trade data
if param.cong==true
    Cjn=results.Djn;
else
    Cjn=results.Cjn;
end

% Construct trade matrix
bilateral_trade_mat=zeros(nb_NUTS,nb_NUTS);
for j=1:g.J
    for n=1:param.N-1
        bilateral_trade_mat( good_to_nuts(n),g.region_id(j) ) = bilateral_trade_mat( good_to_nuts(n),g.region_id(j) )...
            +results.Pjn(j,diff_goods(n) )*Cjn(j,diff_goods(n) );
        
        % Note: I'm using consumer prices, so transport costs are counted in
    end
end

intra_reg_trade_share_differentiated=sum(diag(bilateral_trade_mat))/sum(bilateral_trade_mat(:));


% ----------------
% HOMOGENEOUS GOOD
% ----------------

tol_CV = 1e-4;
MAX_ITER = 30;
n = homogeneous_index; % index of the homogeneous good
intra_reg_trade_share=intra_reg_trade_share_differentiated;

if homogeneous_index~=0 % if there is an homogenenous good
    
    % Toggle proportional vs. domestic priority
    proportional=true;
    
    % Recover trade data
    if param.cong==true
        Cjn=results.Djn;
    else
        Cjn=results.Cjn;
    end
    
    % Compute the diagonal terms (constant)
    Cdiag=zeros(param.J,param.J);
    Tdiag=zeros(param.J,param.J);
    for j=1:param.J
        Tdiag(j,j)=results.Yjn(j,n);
        if proportional==false
            Cdiag(j,j)=min(Cjn(j,n),results.Yjn(j,n));
        else
            total_supply=results.Yjn(j,n);
            for k_id=1:length(g.nodes{j}.neighbors)
                k=g.nodes{j}.neighbors(k_id);
                total_supply=total_supply+results.Qjkn(k,j,n);
            end
            Cdiag(j,j)=results.Yjn(j,n)/total_supply*Cjn(j,n);
        end
    end
    
    % Iterate on system
    C0=Cdiag;
    T0=Tdiag;
    
    
    has_converged=false;
    counter=0;
    while has_converged==false && counter<MAX_ITER
        C1=Cdiag;
        T1=Tdiag;
        
        % Compute nondiagonal terms
        
        for i=1:param.J
            for k=1:param.J
                if i~=k
                    total_imports_k=0;
                    
                    for j_id=1:length(g.nodes{k}.neighbors)
                        j=g.nodes{k}.neighbors(j_id);
                        total_shipments_j=0;
                        for l_id=1:length(g.nodes{j}.neighbors)
                            l=g.nodes{j}.neighbors(l_id);
                            transport_cost=g.delta_tau(j,l).*results.Qjkn(j,l,n)^(1+param.beta)/g.avI(j,l)^param.gamma;
                            total_shipments_j=total_shipments_j+results.Qjkn(j,l,n)+transport_cost;
                        end
                        if total_shipments_j>0
                            T1(i,k)=T1(i,k)+results.Qjkn(j,k,n)/total_shipments_j*(T0(i,j)-C0(i,j));
                        end
                        
                        total_imports_k=total_imports_k+results.Qjkn(j,k,n);
                    end
                    
                    if proportional==false
                        if total_imports_k>0
                            C1(i,k)=T0(i,k)/total_imports_k*(Cjn(k,n)-C0(k,k));
                        end
                    else
                        if total_imports_k+results.Yjn(k,n)>0
                            C1(i,k)=T0(i,k)/(results.Yjn(k,n)+total_imports_k)*Cjn(k,n);
                        end
                    end
                    
                end
                
            end
        end
        
        % Test convergence
        has_converged=max(abs(C1(:)-C0(:))+abs(T1(:)-T0(:)))<tol_CV;
        
        % Update
        C0=C1;
        T0=T1;
        counter=counter+1;
    end
    
    % Store results
    C=C1;
    T=T1;
    
    if ~has_converged
        fprintf('%s.m: could not construct trade matrix.\n',mfilename);
    end
    
    % Aggregate to the NUTS level
    for j=1:param.J
        for k=1:param.J
            bilateral_trade_mat( g.region_id(j), g.region_id(k) ) = bilateral_trade_mat( g.region_id(j), g.region_id(k) )...
                + results.Pjn(k,n)*C(j,k);
            % Note: again, consumer prices
        end
    end
    
    intra_reg_trade_share=sum(diag(bilateral_trade_mat))/sum(bilateral_trade_mat(:));
end
