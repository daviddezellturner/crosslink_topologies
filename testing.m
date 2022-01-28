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

graphs = MMSGraphList;
propogation_times = [];
cycles = cell(1,length(graphs));
bottlenecks = [];
costs_breakdown = cell(1,length(graphs));
cost_calc = [];
lengths = [];

steps = 0:0;%length(graphs)-1;
keys=graphs.keys();
stepsize=keys{2};
steps = steps*stepsize;
for i=steps
    g=distsToRates(graphs(i))
    %[greedy_cycle, cost] = cyclefinder_backtrack_heuristic(g,1);
    [greedy_cycle, cost] = cyclefinder_small(g);
    greedy_costs=[];
    for node = 1:length(greedy_cycle)-1
        greedy_costs = [greedy_costs g(greedy_cycle(node),greedy_cycle(node+1))];
    end
    cycles(1,i/stepsize+1)={greedy_cycle};
    costs_breakdown(1,i/stepsize+1)={greedy_costs};
    bottlenecks = [bottlenecks max(greedy_costs)];
    cost_calc = [cost_calc max(greedy_costs)*length(greedy_cycle)];
    lengths = [lengths length(greedy_cycle)];
end

avg_cost = mean(cost_calc)
avg_length = mean(lengths)
avg_bottleneck = mean(bottlenecks)

%compute travel time as max-cost link times number of links (?)
%explore routing around "problem links"


% Testing small cyclefinder


% plotting
% coordinates

% highlight function
% videowriter