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
    0,0,1,0,0;
    1,0,0,0,0];

%[greedy_cycle, cost] = cyclefinder_backtrack_v2(g_hamiltonian_path_test,3)

% Testing Hamiltonian path algo w/ data propogation on LunaNetWithLLO

lunanetgraphs = dataToGraph("LunaNetWithLLOAccessData.xlsx", 500);
propogation_times = [];
cycles = cell(1,length(lunanetgraphs));

steps = 0:length(lunanetgraphs)-1;
steps = steps*60;
for i=steps
    g=lunanetgraphs(i);
    [greedy_cycle, cost] = cyclefinder_backtrack(g,1);
    t = cyclicPropogationTime(g);
    if length(greedy_cycle) > 1
        propogation_times = [propogation_times t];
    end
    cycles(1,i/60+1)={greedy_cycle};
end

t_avg = mean(propogation_times)