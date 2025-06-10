function flaw_struct = compute_fundamental_law(param,graph,res_baseline,verbose)

if nargin<=3
    verbose=true; % by default, display results
end

% --- Parameters --- %

nperturb=25;

% --- Run replications --- %

logQ_data=zeros(nperturb,1); % local quantity data
logPQ_data=zeros(nperturb,1); % local value data
logI_data=zeros(nperturb,1);
agg_logQ_data=zeros(nperturb,1); % aggregate quantity data
agg_logPQ_data=zeros(nperturb,1); % aggregate value data
agg_logI_data=zeros(nperturb,1);


Qjk_baseline=sum(res_baseline.Qjkn,3); % sum of quantities of flow, no adjustment for value as in Duranton & Tuner
Qjk_baseline=Qjk_baseline+Qjk_baseline';
PQjk_baseline=sum(permute(repmat(res_baseline.Pjn,[1 1 graph.J]),[1 3 2]).*res_baseline.Qjkn,3); % sum value of flows
PQjk_baseline=PQjk_baseline+PQjk_baseline';

if verbose
    fprintf('------------------------------------------------\n');
    fprintf('Evaluating fundamental law of road congestion...\n');
end

for iperturb=1:nperturb
    
    if verbose
        str=sprintf('Computing perturbation %d/%d',iperturb,nperturb);
        fprintf(str);
    end
    
    link=randi([1 graph.ndeg]); % take one link at random
    
    Ijk=res_baseline.Ijk;
    n_to_jk=find(tril(graph.adjacency)==1);
    Ijk(n_to_jk(link))=Ijk(n_to_jk(link))+2*rand();
    Ijk=(Ijk+Ijk')/2; % make it symmetric
    
    [res,flag] = solve_allocation(param,graph,Ijk);
    if ~any( flag.status == [0,1] ) % if IPOPT failed        
        logQ_data(iperturb)=NaN;
        logI_data(iperturb)=NaN;
    else         
        % local quantity data
        Qjk=sum(res.Qjkn,3);
        Qjk=Qjk+Qjk'; % compute total flows in both directions
        
        logQ_data(iperturb)=log(Qjk(n_to_jk(link)))-log(Qjk_baseline(n_to_jk(link)));
        logI_data(iperturb)=log(Ijk(n_to_jk(link)))-log(res_baseline.Ijk(n_to_jk(link)));
        
        % local value data
        PQjk=sum(permute(repmat(res.Pjn,[1 1 graph.J]),[1 3 2]).*res.Qjkn,3);
        PQjk=PQjk+PQjk'; % compute total flows in both directions
                
        logPQ_data(iperturb)=log(PQjk(n_to_jk(link)))-log(PQjk_baseline(n_to_jk(link)));
        
        % aggregate quantity data        
        agg_logQ_data(iperturb)=log(sum(Qjk(:)))-log(sum(Qjk_baseline(:)));
        agg_logI_data(iperturb)=log(sum(Ijk(:)))-log(sum(res_baseline.Ijk(:)));
        
        % aggregate value data                
        agg_logPQ_data(iperturb)=log(sum(PQjk(:)))-log(sum(PQjk_baseline(:)));        
    end
    
    if verbose
        fprintf(repmat('\b',[1 length(str)]));
    end
end


% --- Run the regressions --- %

% Local, quantities

Y=logQ_data(~isnan(logQ_data));
X=logI_data(~isnan(logQ_data));

flaw_struct.local_elasticity_Q_to_I=regress(Y,X);

% Local, values

Y=logPQ_data(~isnan(logQ_data));
flaw_struct.local_elasticity_PQ_to_I=regress(Y,X);

% Aggregate, quantities

Y=agg_logQ_data(~isnan(logQ_data));
X=agg_logI_data(~isnan(logQ_data));

flaw_struct.agg_elasticity_Q_to_I=regress(Y,X);

% Aggregate, values

Y=agg_logPQ_data(~isnan(logQ_data));

flaw_struct.agg_elasticity_PQ_to_I=regress(Y,X);

% --- Report result --- %

if verbose
    fprintf('Local elasticity of logQ to logI = %.4f\n',flaw_struct.local_elasticity_Q_to_I);
    fprintf('Local elasticity of logPQ to logI = %.4f\n',flaw_struct.local_elasticity_PQ_to_I);
    fprintf('Aggregate elasticity of logQ to logI = %.4f\n',flaw_struct.agg_elasticity_Q_to_I);
    fprintf('Aggregate elasticity of logPQ to logI = %.4f\n',flaw_struct.agg_elasticity_PQ_to_I);
end