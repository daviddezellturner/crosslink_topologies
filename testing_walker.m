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

graphs = GPSWalkerApproxGraphList;%MMSGraphList;
propogation_times = [];
cycles = cell(1,length(graphs));
bottlenecks = [];
costs_breakdown = cell(1,length(graphs)-1);
cost_calc = [];
lengths = [];

steps = 1:1;%1:length(graphs)-1;
keys=graphs.keys();
stepsize=keys{2};
steps = steps*stepsize;

cycles = cell(0);
costs = [];
lengths = [];

for i=steps
    
    g=distsToRates(graphs(i));
    [greedy_cycle, cost] = cyclefinder_backtrack_heuristic(g,1)
    g
    %[greedy_cycle, cost] = cyclefinder_small(g)
    cycles{end+1} = greedy_cycle;
    costs = [costs cost];
    lengths = [lengths length(greedy_cycle)];
end

avg_cost = mean(costs);
avg_length = mean(lengths);


% Testing small cyclefinder

v = VideoWriter("walker.avi"); 
v.FrameRate=10;
open(v); 

% plotting
% coordinates

% highlight function
% videowriter

num_nodes=4;

%fps
rate = 10;
t = 0;

even_bottleneck = 0;
odd_bottleneck = 0;

for i=1:floor((length(greedy_cycle)-1)/2)
    if g(greedy_cycle(2*i-1),greedy_cycle(2*i)) > even_bottleneck
        even_bottleneck = g(greedy_cycle(2*i-1),greedy_cycle(2*i));
    end
end

for i=1:floor((length(greedy_cycle)-2)/2)
    if g(greedy_cycle(2*i),greedy_cycle(2*i+1)) > odd_bottleneck
        odd_bottleneck = g(greedy_cycle(2*i),greedy_cycle(2*i+1));
    end
end

t=0;
for i=1:length(greedy_cycle)-1
    i
    if mod(i,2) == 1
        timesteps_here = i*even_bottleneck*rate;
    else
        timesteps_here = i*odd_bottleneck*rate;
    end
    timesteps_here
    for j=1:timesteps_here
        t = t + 1/rate;
        gr = graph(g);
        p = plot(graph(g));
        if mod(i,2) == 1
            for k=1:floor((length(greedy_cycle)-1)/2)
                if k ~= 1
                    highlight(p, [greedy_cycle(2*k-1),greedy_cycle(2*k)], 'EdgeColor','r','LineWidth',1.5);
                end
            end
        else
            for k=1:floor((length(greedy_cycle)-2)/2)
                if j ~= floor((length(greedy_cycle)-2)/2)
                    highlight(p, [greedy_cycle(2*k),greedy_cycle(2*k+1)], 'EdgeColor','r','LineWidth',1.5);
                end
            end
        end
        title("GPS Walker Cyclic Message-Passing: t=" + t + "s");
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        f = getframe(ax,rect);
        writeVideo(v,f);
    end
end

close(v);