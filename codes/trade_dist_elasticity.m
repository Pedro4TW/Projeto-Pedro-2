% compute trade-dist elasticity

function trade_dist_elas = trade_dist_elasticity( results,param,g,GDPj )

        Pmat=repmat(permute(results.Pjn(:,2:end),[3 1 2]),[g.J 1 1]); % use final destination prices
        Q=Pmat.*results.Qjkn(:,:,2:end); % retain only differentiated good
        exports=permute(sum(Q,1),[2 3 1]);
        [~,good_producer_id] = max(param.Zjn(:,2:end),[],1); % provides id of producer per good (inverse of produced_good_id)
        
        Y=[];
        X=[];
        for n=1:param.N-1 % collect data for regression
            set_dest=setdiff(1:g.J,good_producer_id(n))';
            set_dest=set_dest( exports(set_dest,n)>0 ); % retain only places with positive exports
            set_dest=set_dest( g.all_distances(good_producer_id(n),set_dest) > 0 ); % keep distances > 0
            Y=[Y;log(exports(set_dest,n))];
            X=[X;ones(length(set_dest),1),log(GDPj(good_producer_id(n)))*ones(length(set_dest),1),log(GDPj(set_dest)),log(g.all_distances(set_dest,good_producer_id(n)))];
        end
        
        % run regression, compute elasticity
        b=regress(Y,X); % same as X\Y
        trade_dist_elas = b(4);