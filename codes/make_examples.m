%{

    Script to generate the illustrative examples graphs from section 4.
    

%}

clear;
close all;
graph_path = '../draft/graphs_revision/'; % path to save graphs for the paper



%% 4.1.1 One Good on a Regular Geometry
%. COMPARATIVE STATICS OVER K IN A SYMMETRIC NETWORK

% ==============
% INITIALIZATION
% ==============

K=[1,100];
print_dummy=[1,1];
param = init_parameters();
[param,g]=create_graph(param,9,9,'Type','map');

% Define fundamentals
param.N = 1;
param.Zjn = 0.1*ones(g.J,1); % matrix of productivity
Ni=find_node(g,5,5); % center
param.Zjn(Ni)=1; % more productive node

%% Plot the mesh
s='simple_geography_mesh_population';
fig=figure('Units','inches','Position',[0,0,7.5,7.5],'Name',s);
set(gcf,'PaperUnits',get(gcf,'Units'),'PaperPosition',get(gcf,'Position')+[-.25 2.5 0 0]);
plot_graph(param,g,[],'Edges','off','Mesh','on','MeshColor',[.2 .4 .6]);
box off;
axis off;
print('-depsc2',[graph_path,s,'.eps']);
saveas(fig,[graph_path,s,'.jpg']);

s='simple_geography_mesh_productivity';
fig=figure('Units','inches','Position',[0,0,7.5,7.5],'Name',s);
set(gcf,'PaperUnits',get(gcf,'Units'),'PaperPosition',get(gcf,'Position')+[-.25 2.5 0 0]);
sizes=ones(g.J,1);
sizes(find_node(g,5,5))=2.5;
plot_graph(param,g,[],'Edges','off','Shades',(param.Zjn-.1)/.9,'Sizes',sizes,'NodeFgColor',[.6 .8 1],'Mesh','on','MeshColor',[.2 .4 .6]);
box off;
axis off;
print('-depsc2',[graph_path,s,'.eps']);
saveas(fig,[graph_path,s,'.jpg']);

%% Compute networks
for i=1:length(K) % Solve for the optimal network for each K
    param.K=K(i);
    results(i,1)=optimal_network(param,g);
end

%% Plot networks

for k=1:length(K)
    if print_dummy(k)==true
        s=sprintf('Comparative_statics_K=%d',K(k));

        fig=figure('Units','inches','Position',[0,0,7.5,6],'Name',s);
        set(gcf,'PaperUnits',get(gcf,'Units'),'PaperPosition',get(gcf,'Position')+[-.25 2.5 0 0]);
        set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times');

        hor_scaling=1.1;
        ver_scaling=1.05;

        h=subplot(2,2,1);
        size=.25*ones(g.J,1);
        plot_graph( param,g,results(k,1).Ijk,'Sizes',size );
        pos=get(h,'Position');
        set(h,'Position',[pos(1) pos(2) hor_scaling*pos(3) ver_scaling*pos(4)]);
        title('(a) Transport Network (I_{jk})','FontWeight','normal','FontName','Times','FontSize',9);

        h=subplot(2,2,2);
        plot_graph( param,g,results(k,1).Qjkn,'Nodes','off','Arrows','on','ArrowScale',1 );
        pos=get(h,'Position');
        set(h,'Position',[pos(1) pos(2) hor_scaling*pos(3) ver_scaling*pos(4)]);
        title('(b) Shipping (Q_{jk})','FontWeight','normal','FontName','Times','FontSize',9);

        h=subplot(2,2,3);
        plot_graph( param,g,results(k,1).Ijk,'Nodes','off','Map', results(k,1).Pjn/max(results(k,1).Pjn) );
        pos=get(h,'Position');
        set(h,'Position',[pos(1) pos(2) hor_scaling*pos(3) ver_scaling*pos(4)]);
        colorbar();
        title('(c) Prices (P_{j})','FontWeight','normal','FontName','Times','FontSize',9);

        h=subplot(2,2,4);
        plot_graph( param,g,results(k,1).Ijk,'Nodes','off','Map',results(k,1).cj/max(results(k,1).cj) );
        hold off;
        pos=get(h,'Position');
        set(h,'Position',[pos(1) pos(2) hor_scaling*pos(3) ver_scaling*pos(4)]);
        colorbar();
        title('(d) Consumption (c_{j})','FontWeight','normal','FontName','Times','FontSize',9);

        print('-depsc2',[graph_path,s,'.eps'],'-r400');
        saveas(fig,[graph_path,s,'.jpg']);
    end
end


% 4.1.2 RANDOM CITIES
%
% Init and Solve network

width=9; height=9;
nb_cities=20;

param = init_parameters('TolKappa',1e-5,'K',100,'Annealing','off');
[param,g]=create_graph(param,width,height,'Type','triangle'); % case with random cities, one good

% set fundamentals

rng(5);
param.N = 1;
param.Zjn = param.minpop*ones(g.J,1); % matrix of productivity
param.Hj = ones(g.J,1); % matrix of housing
param.Lj = 0*ones(g.J,1); % matrix of population

Ni=find_node(g,ceil(width/2),ceil(height/2)); % find center
param.Zjn(Ni)=1; % more productive node
param.Lj(Ni)=1; % more productive node

for i=1:nb_cities-1
    newdraw=false;
    while newdraw==false
        j=nearest(1+rand()*(g.J-1));
        if param.Lj(j)<=param.minpop
            newdraw=true;
            param.Lj(j)=1;
        end
    end
end

param.hj=param.Hj./param.Lj;
param.hj(param.Lj==0)=1;

% Convex case
results(1)=optimal_network(param,g);

% Nonconvex - no annealing
param.gamma=2;
results(2)=optimal_network(param,g);

% Nonconvex - annealing
results(3)=annealing(param,g,results(2).Ijk,'PerturbationMethod','rebranching');

welfare_increase=consumption_equivalent(param,g,results(2).cj,results(2).Lj,results(3).welfare)

% plot
close all;

s={'random_cities_convex','random_cities_nonconvex','random_cities_annealing'};
titles = [ {'Transport Network (I_{jk})','Shipping (Q_{jk})'};
           {'Transport Network (I_{jk})','Shipping (Q_{jk})'};
           {'Before annealing','After annealing'} ];

plots = [ {'results(i).Ijk','results(i).Qjkn'};
          {'results(i).Ijk','results(i).Qjkn'};
          {'results(2).Ijk','results(3).Ijk'}];

texts = [ {'',''};
         {'',''};
         {'Welfare = 1',sprintf('Welfare = %1.3f (cons. eq.)',welfare_increase)}];

arrows = [ {'off','on'};
           {'off','on'};
           {'off','off'} ];


for i=1:3
    fig=figure('Units','inches','Position',[0,0,7.5,3],'Name',char(s(i)));

    subplot(1,2,1);
    plot_graph(param,g,eval(char(plots(i,1))),'Sizes',1.2*param.Lj,'Arrows',char(arrows(i,1)));
    title(sprintf('(%c) %s',97+2*(i-1),char(titles(i,1))),'FontWeight','normal','FontName','Times','FontSize',9);
    text(.5,-.05,char(texts(i,1)),'HorizontalAlignment','center','Fontsize',8);


    subplot(1,2,2);
    plot_graph(param,g,eval(char(plots(i,2))),'Sizes',1.2*param.Lj,'Arrows',char(arrows(i,2)),'ArrowScale',1);
    title(sprintf('(%c) %s',97+2*(i-1)+1,char(titles(i,2))),'FontWeight','normal','FontName','Times','FontSize',9);
    text(.5,-.05,char(texts(i,2)),'HorizontalAlignment','center','Fontsize',8);


    print('-depsc2',[graph_path,char(s(i)),'.eps'],'-r400');
    saveas(fig,[graph_path,char(s(i)),'.jpg']);

end

%% ------------------------------------
% 4.2 Random Cities with Multiple Goods
% -------------------------------------

clear results;

ngoods=10;
width=9; height=9;

% Init and Solve network

param = init_parameters('TolKappa',1e-5,'K',10000,'LaborMobility','on','Annealing','off');
[param,g]=create_graph(param,width,height,'Type','triangle');

% set fundamentals

param.N = 1+ngoods;
param.Zjn = 0*ones(g.J,param.N); % matrix of productivity

param.Zjn(:,1)=1; % the entire countryside produces the homogenenous good 1

if ngoods>0
    Ni=find_node(g,5,5); % center
    param.Zjn(Ni,2)=1; % central node
    param.Zjn(Ni,1)=0; % central node
    
    rng(5);
    for i=2:ngoods
        newdraw=false;
        while newdraw==false
            j=round(1+rand()*(g.J-1));
            if param.Zjn(j,1)>0
                newdraw=true;
                param.Zjn(j,i+1)=1;
                param.Zjn(j,1)=0;
            end
        end
    end
end

% Convex case
results(1)=optimal_network(param,g);

% Nonconvex
param.gamma=2;
results(2)=optimal_network(param,g,results(1).Ijk);


%% Plot results
close all;

cols = 3; %4 % number of columns
rows = ceil((1+param.N)/cols);

s={'random_cities_multigoods_convex','random_cities_multigoods_nonconvex'};

for j=1:2
    fig=figure('Units','inches','Position',[0,0,7.5,11],'Name',char(s(j)));
%     fig=figure('Units','inches','Position',[0,0,11,7.5],'Name',char(s(j)));  % for horizontal display


    % Plot network
    subplot(rows,cols,1);
    plot_graph( param,g,results(j).Ijk, 'Shades', (results(j).Lj-min(results(j).Lj))./(max(results(j).Lj)-min(results(j).Lj)), 'Sizes', 1+16*(results(j).Lj./mean(results(j).Lj)-1), 'NodeFgColor',[.6 .8 1],'Transparency','off' );
    h=title('(a) Transport Network (I_{jk})','FontWeight','normal','FontName','Times','FontSize',9);

    for i=1:param.N
        subplot(rows,cols,i+1);
        sizes=3*results(j).Yjn(:,i)./sum(results(j).Yjn(:,i));
        shades=param.Zjn(:,i)./max(param.Zjn(:,i));
        plot_graph( param,g,results(j).Qjkn(:,:,i),'Arrows','on','ArrowScale',1,'ArrowStyle','thin','Nodes','on','Sizes',sizes,'Shades',shades, 'NodeFgColor',[.6 .8 1],'Transparency','off' );
        title(sprintf('(%c) Shipping (Q^{%d}_{jk})',97+i,i),'FontWeight','normal','FontName','Times','FontSize',9);
    end

    % Save
    print('-depsc2',[graph_path,char(s(j)),'.eps'],'-r400');
    saveas(fig,[graph_path,char(s(j)),'.jpg']);
end

%% ------------------------------------
% 4.3.1 Geography
% -------------------------------------

clear results;

width=13; height=13;

% Init and Solve network

param = init_parameters('TolKappa',1e-5,'K',100,'LaborMobility','off');
[param,g]=create_graph(param,width,height,'Type','map'); 

% set fundamentals

rng(5);
param.N = 1;
param.Zjn = param.minpop*ones(g.J,1); % matrix of productivity
param.Hj = ones(g.J,1); % matrix of housing
param.Lj = 0*ones(g.J,1); % matrix of population

Ni=find_node(g,5,5); % center
param.Zjn(Ni)=1; % more productive node
param.Lj(Ni)=1; % more productive node

nb_cities=20; % draw a number of random cities in space
for i=1:nb_cities-1
    newdraw=false;
    while newdraw==false
        j=round(1+rand()*(g.J-1));
        if param.Lj(j)<=param.minpop
            newdraw=true;
            param.Lj(j)=1;
        end
    end
end

param.hj=param.Hj./param.Lj;
param.hj(param.Lj==0)=1;

% --------------
% Draw geography

z=zeros(g.J,1); % altitude of each node
obstacles = [];

geography = struct('z',z,'obstacles',obstacles);
g = apply_geography( g, geography );

param0=param; % store initial params
g0=g; % store initial graph


%% Blank geography
results(1)=optimal_network(param,g);

%% Mountain
mountain_size=.75;
mountain_height=1;
mount_x=10;
mount_y=10;
geography(2)=geography(1);
geography(2).z=mountain_height*exp(-((g.x-mount_x).^2 + (g.y-mount_y).^2)/(2*mountain_size^2));

g = apply_geography( g, geography(2),'AlphaUp_i',10,'AlphaDown_i',10 );
results(2)=optimal_network(param,g);

%% Adding river and access by land
g=g0;
geography(3)=geography(2);
geography(3).obstacles = [6 + (1-1)*width, 6 + (2-1)*width;
    6 + (2-1)*width, 6 + (3-1)*width;
    6 + (3-1)*width, 7 + (4-1)*width;
    7 + (4-1)*width, 8 + (5-1)*width;
    8 + (5-1)*width, 9 + (5-1)*width;
    11 + (5-1)*width, 12 + (5-1)*width;
    12 + (5-1)*width, 13 + (5-1)*width]; % Nobj x 2 matrix of (i,j) pairs of locations where a geographic barrier should be drawn

g = apply_geography( g, geography(3),'AcrossObstacleDelta_i', inf, 'AlphaUp_i',10,'AlphaDown_i',10 );
results(3)=optimal_network(param,g);

%% Reinit and put another river and bridges
g=g0;
geography(4)=geography(1);
geography(4).z=mountain_height*exp(-((g.x-mount_x).^2 + (g.y-mount_y).^2)/(2*mountain_size^2));

geography(4).obstacles = [6 + (1-1)*width, 6 + (2-1)*width;
    6 + (2-1)*width, 6 + (3-1)*width;
    6 + (3-1)*width, 7 + (4-1)*width;
    7 + (4-1)*width, 8 + (5-1)*width;
    8 + (5-1)*width, 9 + (5-1)*width;
    9 + (5-1)*width, 10 + (5-1)*width;
    10 + (5-1)*width, 11 + (5-1)*width;
    11 + (5-1)*width, 12 + (5-1)*width;
    12 + (5-1)*width, 13 + (5-1)*width];

g = apply_geography( g, geography(4), 'AlphaUp_i',10,'AlphaDown_i',10, 'AcrossObstacleDelta_i',2,'AlongObstacleDelta_i',inf );
results(4)=optimal_network(param,g);

%% Allowing for water transport
g=g0;
geography(5)=geography(1);
geography(5).z=mountain_height*exp(-((g.x-mount_x).^2 + (g.y-mount_y).^2)/(2*mountain_size^2));
geography(5).obstacles = [6 + (1-1)*width, 6 + (2-1)*width;
    6 + (2-1)*width, 6 + (3-1)*width;
    6 + (3-1)*width, 7 + (4-1)*width;
    7 + (4-1)*width, 8 + (5-1)*width;
    8 + (5-1)*width, 9 + (5-1)*width;
    9 + (5-1)*width, 10 + (5-1)*width;
    10 + (5-1)*width, 11 + (5-1)*width;
    11 + (5-1)*width, 12 + (5-1)*width;
    12 + (5-1)*width, 13 + (5-1)*width];

g = apply_geography( g, geography(5), 'AlphaUp_i',10,'AlphaDown_i',10, 'AcrossObstacleDelta_i',2,'AlongObstacleDelta_i',.5 );
results(5)=optimal_network(param,g);

%% Increasing returns to transport
param.gamma=2;
geography(6)=geography(5);
results(6)=optimal_network(param,g);

%% Plot results
close all;

s={'geography_blank','geography_mountain','geography_river','geography_bridges','geography_water_transport','geography_increasing_returns'};
obstacles = {'off','off','on','on','on','on'};

for j=1:length(s)
    fig=figure('Units','inches','Position',[0,0,7.5,5],'Name',char(s(j)));

    % Plot network
    plot_graph( param,g,results(j).Ijk,'Geography','on','GeographyStruct',geography(j),'Sizes',3*results(j).cj.*(param.Lj>param.minpop)/max(results(j).cj),'Mesh','on','MeshTransparency',.2,'Obstacles',char(obstacles(j)), 'Shades',param.Zjn, 'MinEdgeThickness', 1.5 );

    % Save
    print('-depsc2',[graph_path,char(s(j)),'.eps']);
    saveas(fig,[graph_path,char(s(j)),'.png']);
end
