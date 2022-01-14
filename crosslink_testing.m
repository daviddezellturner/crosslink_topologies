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
%_cycle, cost] = cyclefinder(g,3);


%% TEST PATHFINDER
%[best_path, cost] = pathfinder(g, 3, 1);

%% TEST GENETIC
% parameters = genetic(frames)

%% TEST METRICS
% (figures of merit)

%% TEST LINKSTATE
linkstate(g)
