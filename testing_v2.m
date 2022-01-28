%% GOALS
%% 1. Test efficacy of greedy path vs tree propogation vs ... for various network types
%% 2. Test/tune paramaters of individual path-finding algos to optimize for metrics
%% 3. Do this for orbit propogation data, not just toy examples as below

% graph g consists of three nodes.
% links:
% 1 -> 2 with weight 3
% 1 -> 3 with weight 7
% 2 -> 1 with weight 1
% 2 -> 3 with weight 2
% 3 -> 1 with weight 9]
% no link from 3 -> 2
% no self-links (we input weight zero if there is no link.)
g = [0, 3, 7; 1, 0, 2; 9, 0, 0];

% toy example for genetic algo: consider 2 timesteps of a 6-satellite, 2-ring structure
% pretend each config is attained around half of the time (see sketches)
% initialize simulated distance-finding costs as a function of actual costs
% and tune function parameters to give best performance according to all metrics
G1 = [0,0,1,1,0,0;
    0,0,0,0,1,1;
    1,0,0,1,0,1;
    1,0,1,0,1,0;
    0,1,0,1,0,1;
    0,1,1,0,1,0];
G2 = [0,0,1,0,0,0;
    0,0,0,0,1,0;
    1,0,0,1,0,1;
    0,0,1,0,1,0;
    0,1,0,1,0,1;
    0,0,1,0,1,0];
frames = [G1, G2];

%% TEST CYCLEFINDER
g_backtrack_test = [0,1,0,0,0;
    0,0,1,10,0;
    0,0,0,1,1;
    0,0,1,0,0;
    1,0,0,0,0];

g_hamiltonian_path_test = [0,1,0,0,0;
    0,0,1,10,0;
    0,0,0,1,1;
    0,0,1,08,0;
    1,0,0,0,0];

%[greedy_cycle, cost] = cyclefinder_backtrack_v2(g_hamiltonian_path_test,3)

% Testing Hamiltonian path algo w/ data propogation on LunaNetWithLLO

graphs = LunaNetWithLLOGraphList;
propogation_times = [];
cycles = cell(1,length(graphs));
bottlenecks = [];
costs_breakdown = cell(1,length(graphs));
cost_calc = [];
lengths = [];

steps = 0:200;
keys=graphs.keys();
stepsize=keys{2};
steps = steps*stepsize;

cycles = cell(0);
costs = [];
lengths = [];

for i=steps
    
    g=distsToRates(graphs(i));
    %[greedy_cycle, cost] = cyclefinder_backtrack_heuristic(g,1)
    [greedy_cycle, cost] = cyclefinder_small(g)
    cycles{end+1} = greedy_cycle;
    costs = [costs cost];
    lengths = [lengths length(greedy_cycle)];
end

avg_cost = mean(costs);
avg_length = mean(lengths);


% Testing small cyclefinder
% 
% v = VideoWriter("mms.avi"); 
% v.FrameRate=10;
% open(v); 
% 
% % plotting
% % coordinates
% 
% % highlight function
% % videowriter
% 
% num_nodes=4;
% 
% %fps
% rate = 10;
% t = 0;
% for i=1:length(greedy_cycle)-1
%     timesteps_here = round(min(i,num_nodes)*g(greedy_cycle(i),greedy_cycle(i+1))*rate);
%     for j=1:timesteps_here
%         t = t + 1/rate;
%         gr = graph(g);
%         p = plot(graph(g),'EdgeLabel',graph(g).Edges.Weight);
%         highlight(p, [greedy_cycle(i),greedy_cycle(i+1)], 'EdgeColor','r','LineWidth',1.5);
%         title("MMS Cyclic Message-Passing: t=" + t + "s");
%         ax = gca;
%         ax.Units = 'pixels';
%         pos = ax.Position;
%         ti = ax.TightInset;
%         rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
%         f = getframe(ax,rect);
%         writeVideo(v,f);
%     end
% end
% for i=1:length(greedy_cycle)-1
%     num_nodes*g(greedy_cycle(i),greedy_cycle(i+1))
%     timesteps_here = round(num_nodes*g(greedy_cycle(i),greedy_cycle(i+1))*rate);
%     num_nodes*g(greedy_cycle(i),greedy_cycle(i+1))
%     for j=1:timesteps_here
%         t = t + 1/rate;
%         gr = graph(g);
%         p = plot(graph(g),'EdgeLabel',graph(g).Edges.Weight);
%         highlight(p, [greedy_cycle(i),greedy_cycle(i+1)], 'EdgeColor','r','LineWidth',1.5);
%         title("MMS Cyclic Message-Passing: t=" + t + "s");
%         ax = gca;
%         ax.Units = 'pixels';
%         pos = ax.Position;
%         ti = ax.TightInset;
%         rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
%         f = getframe(ax,rect);
%         writeVideo(v,f);
%     end
% end
% 
% close(v);