% calibrate.m: this script calibrates:
% if 1) mobile labor: productivity parameter Zjn and housing Hj
%       in order to match GDP and population for all locations.
%    2) fixed labor: productivity parameter Zjn in order to match GDP for all locations.
%    3) if calibrate_d0 == true, calibrates d0 to match intraregional trade
%
% WARNING: this only works with exactly one good produced by location.

% -----------------------------
% Parameters of the calibration

has_calibration_converged=false;
counter=0;
tol_L=1e-2;
tol_GDP_per_cap=1e-2;
tol_intra_reg_trade_share=0.5e-2;
tol_trade_dist_elas=1e-2; % put 3e-2 if target elasticity -1
MAX_ITER=400;
step_Z=0.95;
step_L=0.75;
step_d0=0.95;
step_d1=0.95;
display_details = false;
save_before_it_crashes = false;

if ~exist('calibrate_d0')
    calibrate_d0=false;
end

if ~exist('calibrate_d1')
    calibrate_d1=false;
end

if ~exist('compute_trade_dist_elas')
    compute_trade_dist_elas=true;
end

if ~exist('verbose')
    verbose = false;
end

% --------------
% Define targets

target_L = g.L;
target_GDP_per_cap = g.Y./g.L;
target_intra_reg_trade_share = [0.53, 0.44];% NUTS1: 0.53, NUTS2: 0.44% old: 0.39;
target_intra_reg_trade_share = target_intra_reg_trade_share(NUTS);
target_trade_dist_elas = -1.05; % -1

exclude_locations = target_GDP_per_cap / mean (target_GDP_per_cap) > 2.5; % exclude some locations
avg_dist= mean( g.distance ( g.adjacency == 1) );

% --------------------------------------------
% Identify which good is produced by each city
[~,produced_good_id] = max(param.Zjn,[],2);

% -------------
% Define bounds

% the following are just rough guesses to get better starting point
Cn=sum(param.Zjn.*g.L(:,ones(param.N,1)).^param.a,1); % total quantity produced of each good
Pn=((Cn/Cn(1)).^(-1/param.sigma))';

Zu = 1e2*g.Y./(Pn(produced_good_id).*g.L.^(param.a)); % upper bound on productivity
Zl = 1e-2*g.Y./(Pn(produced_good_id).*g.L.^(param.a)); % lower bound

hu=1e2*ones(g.J,1); % upper bound on housing per capita
hl=1e-2*ones(g.J,1); % lower bound

% define some random starting points
if calibrate_d0==true
    switch param.mobility
        case true
            if beta>0.5
                d0=1;
            else
                d0=0.001;
            end
        case false
            if beta>0.5
                d0=1;
            else
                d0=0.001;
            end
    end
end
if calibrate_d1==true
    switch param.mobility
        case true
            d1=.85;%.85;
        case false
            d1=1;
    end
end

d0u=d0*avg_dist^d1*100; % rescale the d0 by average distance for greater stability
d0l=d0*avg_dist^d1/100;

d1u=d1*1.2; % upper bound on delta1
d1l=d1/1.2; % lower bound


% ------------------------------------------------------
% Run iterations (bisection algorithm on the Zjn and Hj)
% ------------------------------------------------------

x0=[]; % we init optimization from previous optimum, empty for the first run
best_fit.score=inf;


skip_update=false; % dummy indicating whether to skip updating once (debugging)
attempts=0;
Z_MEAN=10; % normalization (nominal GDP is normalized afterwards), chose>1 for numerical reasons
while ~has_calibration_converged && counter<MAX_ITER
    
    % Evaluate at current middle point
    
    Zcurrent=exp(0.5*log(Zu)+0.5*log(Zl));
    Zcurrent=Z_MEAN*Zcurrent/mean(Zcurrent);
    param.Zjn((produced_good_id-1)*g.J+(1:g.J)')=Zcurrent;
    
    hcurrent = exp(0.5*log(hu)+0.5*log(hl));
    param.hj = hcurrent;
    param.Hj = param.hj .* target_L;
    
    d0current=exp(0.5*log(d0u)+0.5*log(d0l));
    d1current=exp(0.5*log(d1u)+0.5*log(d1l));
    
    % Compute tau or kappa given d0current
    g.delta_tau = d0current* ( g.distance ).^d1current / avg_dist^d1current;
    
    % Solve allocation
    if save_before_it_crashes==true
        save(debug_file_str,'param','g','x0','Zu','Zl','Zcurrent','hu','hl','hcurrent','d0u','d0l','d0current','d1u','d1l','d1current');
    end
    
    t0=clock();
    [results,flag,x1] = solve_allocation(param,g,g.avI,verbose,x0);
    t1=clock();
    
    % The following deals with some failure of replicability that we have
    % found using IPOPT:
    if ~any( flag.status == [0,1] )
        
        x1=[]; % try again with default seed
        skip_update=true;
        attempts=attempts+1;
        if attempts==3 % third time it fails, return error
            error('%s.m: IPOPT returned with error code %d.',mfilename,flag.status);
        end
    else
        attempts=0;
    end
    
    if skip_update==true && any( flag.status == [0,1] )
        skip_update=false;
    end
    
    % Compute statistics
    PHj = param.alpha / (1-param.alpha) * results.PCj .* results.hj ./ results.cj;
    GDPj = sum(results.Pjn.*results.Yjn,2)+PHj.*param.Hj;
    GDPj = GDPj/sum(GDPj); % normalize nominal GDP by 1
    GDP_per_cap_j = GDPj./results.Lj;
    GDP_per_cap_j(results.Lj==0) = 0;
    
    
    if isfield(g,'nuts')
        [~,intra_reg_trade_share] = compute_trade_shares(param,g,results,1); % first is all, second is differentiated
    end
    
    % compute trade-distance elasticity
    if compute_trade_dist_elas==true
        trade_dist_elas = trade_dist_elasticity( results,param,g,GDPj );
    else
        trade_dist_elas = 0; % just set trade_dist_elas to 0
    end
    
    % Evaluate fit
    diff_GDP_per_cap=(GDP_per_cap_j-target_GDP_per_cap)./target_GDP_per_cap;
    diff_GDP_per_cap(target_GDP_per_cap==0)=0; % the places with 0 income maintain 0 income throughout
    distance_GDP_per_cap=max(abs(diff_GDP_per_cap(~exclude_locations)));
    score=distance_GDP_per_cap;
    has_GDP_converged=distance_GDP_per_cap<tol_GDP_per_cap;
    has_calibration_converged = has_GDP_converged;    
    
    
    if param.mobility==true
        diff_L=results.Lj-target_L;
        distance_L=max(abs(diff_L(~exclude_locations)));
        score=score+distance_L;
        has_L_converged=distance_L<tol_L;
        
        has_calibration_converged=has_calibration_converged & has_L_converged;
    end
    
    if isfield(g,'nuts')
        diff_intra_reg_trade_share = intra_reg_trade_share-target_intra_reg_trade_share;
        distance_intra_reg_trade_share = abs(diff_intra_reg_trade_share);
        score=score+distance_intra_reg_trade_share;
    end
    if calibrate_d0==true
        has_intra_reg_trade_share_converged = distance_intra_reg_trade_share < tol_intra_reg_trade_share;
        has_calibration_converged=has_calibration_converged & has_intra_reg_trade_share_converged;
    end
    
    if calibrate_d1==true
        diff_trade_dist_elas = trade_dist_elas-target_trade_dist_elas;
        distance_trade_dist_elas = abs(diff_trade_dist_elas);
        score=score+distance_trade_dist_elas;
        has_trade_dist_elas_converged = distance_trade_dist_elas < tol_trade_dist_elas;
        
        has_calibration_converged=has_calibration_converged & has_trade_dist_elas_converged;
    end
    
    % Update bounds
    if ~has_GDP_converged && skip_update==false
        Ipos=find(diff_GDP_per_cap>tol_GDP_per_cap & ~exclude_locations);
        Ineg=find(diff_GDP_per_cap<-tol_GDP_per_cap & ~exclude_locations);
        
        Zu(Ipos)=Zcurrent(Ipos)+step_Z*(Zu(Ipos)-Zcurrent(Ipos)); % bring the bounds closer together slowly
        Zl(Ineg)=Zcurrent(Ineg)-step_Z*(Zcurrent(Ineg)-Zl(Ineg));
    end
    
    if param.mobility==true
        if ~has_L_converged && skip_update==false
            Ipos=find(diff_L>tol_L & ~exclude_locations);
            Ineg=find(diff_L<-tol_L & ~exclude_locations);
            
            hu(Ipos)=hcurrent(Ipos)+step_L*(hu(Ipos)-hcurrent(Ipos)); % bring the bounds closer together slowly
            hl(Ineg)=hcurrent(Ineg)-step_L*(hcurrent(Ineg)-hl(Ineg));
        end
    end
    
    if calibrate_d0==true && skip_update==false
        if ~has_intra_reg_trade_share_converged && distance_GDP_per_cap<0.05
            if diff_intra_reg_trade_share>tol_intra_reg_trade_share
                d0u = d0current + step_d0*(d0u-d0current);
            end
            
            if diff_intra_reg_trade_share<-tol_intra_reg_trade_share
                d0l = d0current - step_d0*(d0current-d0l);
            end
        end
    end
    
    if calibrate_d1==true && skip_update==false
        if ~has_trade_dist_elas_converged && max(abs(diff_GDP_per_cap))<0.05 && distance_intra_reg_trade_share < 0.02
            if diff_trade_dist_elas>tol_trade_dist_elas
                d1l = d1current - step_d1*(d1current-d1l);
            end
            
            if diff_trade_dist_elas<-tol_trade_dist_elas
                d1u = d1current + step_d1*(d1u-d1current);
            end
        end
    end
    
    if score<best_fit.score
        best_fit.results=results;
        best_fit.score=score;
        best_fit.Zjn=param.Zjn;
        best_fit.hj=param.hj;
        best_fit.Hj=param.Hj;
        best_fit.delta_tau=g.delta_tau;
        best_fit.d0=d0current/avg_dist^d1current;
        best_fit.d1=d1current;
        if isfield(g,'nuts')
            best_fit.intra_reg_trade_share = intra_reg_trade_share;
        end
        best_fit.trade_dist_elas = trade_dist_elas;
    end
    
    counter=counter+1;
    x0=x1;
    
    old_score=score;
    old_param=param;
    old_results=results;
    old_GDP_per_cap_j = GDP_per_cap_j;
    
    % Display iteration results
    
    fprintf('CALIBRATION: iteration No.%d...\n',counter);
    vars = {'GDP_per_cap'};
    if param.mobility==true
        vars = [vars;'L'];
    end
    
    if isfield(g,'nuts')
        vars = [vars;'intra_reg_trade_share'];
    end
    
    if calibrate_d1==true
        vars = [vars;'trade_dist_elas'];
    end
    
    str = 'Fit max distance [';
    for k=1:length(vars)
        str = [str,vars{k}];
        if k<length(vars)
            str = [str,' '];
        end
    end
    str = [str,'] = ['];
    for k=1:length(vars)
        str = [str,sprintf('%2.3f',eval(['distance_',vars{k}]))];
        if k<length(vars)
            str = [str,' '];
        end
    end
    str = [str,']'];
    
    disp(str);
    
    if isfield(g,'nuts')
        fprintf('Moments: intra_reg_trade_share=%.3f, trade_dist_elas=%.3f\n',intra_reg_trade_share,trade_dist_elas);
    else
        fprintf('Moments: trade_dist_elas=%.3f\n',trade_dist_elas);
    end
    fprintf('Current estimates: d0=%.5f, d1=%.3f\n',d0current/avg_dist^d1current,d1current);
    if display_details==true
        fprintf('GDP_per_capita fit [model data diff]:\n');
        [GDP_per_cap_j target_GDP_per_cap diff_GDP_per_cap]
        
        if param.mobility==true
            fprintf('Population fit:\n');
            [results.Lj target_L diff_L]
        end
    end
    fprintf('Computation time: %3.1f secs.\n',etime(t1,t0));
end

% store the estimated parameters d0
d0=d0current/avg_dist^d1current;
d1=d1current;

if ~has_calibration_converged % if hasn't converged return best fit
    
    d0=best_fit.d0;
    d1=best_fit.d1;
    if isfield(g,'nuts')
        intra_reg_trade_share = best_fit.intra_reg_trade_share;
    end
    trade_dist_elas = best_fit.trade_dist_elas;
    g.delta_tau = best_fit.delta_tau;
    param.Zjn = best_fit.Zjn;
    param.hj = best_fit.hj;
    param.Hj = best_fit.Hj;
    results = best_fit.results;
    results.Ijk = g.avI;
    
end

% store parameters and targets
param.d0 = d0;
param.d1 = d1;
if calibrate_d0
    results.intra_reg_trade_share = best_fit.intra_reg_trade_share;
end
results.trade_dist_elas = best_fit.trade_dist_elas;
