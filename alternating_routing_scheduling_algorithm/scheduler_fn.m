function [node_schedules, msg_arrivals, Ga_n2a, action_schedule] = scheduler_fn(m_paths, travel_time, plot_on)
% Maryann Rui
% January 2022 

% INPUTS
% m_paths = container of message paths. m_paths(i) = digraph of path/tree 
% that message i takes. 
% travel_time = function(j, k, i) returning time it takes message of size
% m_i kb [incorporated in the function] on edge j->k (may be based on slant
% range distances, also incorporated)
% plot_on = (default false) plot action graph & intermediate fringe graphs; 

% OUTPUTS
% node_schedules = cell array of node schedules: node_schedules{i} is a kx3
% matrix where k = # actions node i performs. row1 = start time, row2 = end
% time, row3 = message being passed 
% msg_arrivals = nxn matrix, msg_arrivals(i, j) = first time at which
% message i arrives at node j according to schedule. NaN if it never does.
% Ga_n2a = n_a x 4 array. each row corresponds to action. col2 = action
% start node, col3 = action end node, col4 = message being passed.
% action_schedule = (# actions)x5 array. each row corresponds 
% to an action executed:[t_start t_end node_start node_end message]
if nargin == 2
    plot_on = false; 
end

%% Outline of the program: 
% 1. Given message paths
% 2. Create action tree with MPOC edges
% 3. Iterate until no more unvisited actions remain: 
%   - Create current fringe graph with NTOC edges and heuristic weights.
%   - Find or estimate the maximum weight independent set (WIS).
%   - Mark actions in this WIS as visited in the main action graph G_a
% 4. Schedule the actions globally. Each action involves two (original
% graph) nodes-- find the earliest (left aligned) time in which these two
% nodes are both free, but after the nodes have executed their previously
% assigned actions (do not change the order of actions even if there is
% room in the schedule. This is due to MPOC contraints). Keep track of the
% time at which each message arrives at its destination node; especially
% the time of the last message to arrive, which is the total message
% propagation time. 
% (5). Estimate message propagation times along given nodes from the
% schedule. \theta^i_jk = \tau^i_k - \tau^i_j [incorporates Queue delays]
%% 1. Given message paths
% > Message paths are given in a dictionary (container): 
% {i: G_i} where i is originating node (sending message i), and G_i is a
% directed graph with one source node (i) and one sink (destination)

%% 2. Create action tree with MPOC edges

% Maintain a map between Ga's node labels (1, 2, ...) and the actions they represent. 
Ga_n2a = []; % n2a = nodes to actions 
% col 1 = Ga node #, cols 2-3 are the edge in G the action corresponds to, 
% and col 4 is the message # sent in this action
Ga_edges_source = []; 
Ga_edges_sink = []; 

nodeind = 1; % keep track of action node index
for mi = keys(m_paths) % for each message path, add MPOC edges in Ga
    mi = mi{1}; % extract node # from cell 
    Gi = m_paths(mi); % get message path graph
    actionsi = Gi.Edges.EndNodes; % get actions (edges in path) kx2
    num_actionsi = size(actionsi, 1);
    % label these actions with their node #s in Ga
    nodenums = nodeind:(nodeind+num_actionsi-1); 
    nodeind = nodeind+num_actionsi; 
    Ga_n2a = [Ga_n2a; nodenums', actionsi , ...
        ones(num_actionsi, 1)*mi]; 
    % MPOC edges in action graph 
    for ind = 1:num_actionsi
        endnode = actionsi(ind, 2); % end node of action. see if it connects with any other actions, which would warrant a MPOC edge
        connected_actions = find(actionsi(:, 1) == endnode); 
        if ~isempty(connected_actions)
            Ga_edges_source = [Ga_edges_source nodenums(ind)*ones(1, numel(connected_actions))]; 
            Ga_edges_sink = [Ga_edges_sink nodenums(connected_actions)]; 
        end
    end
end
num_actions = nodeind-1; 
% create action graph
Ga_adj = sparse(Ga_edges_source, Ga_edges_sink, 1, num_actions, num_actions); 
Ga = digraph(Ga_adj); 
if plot_on
    figure(); 
    plot(Ga); 
    title('action graph G_a'); 
end

%% 3. Iterate: 
%% 3a. Create fringe graph, find max IS, mark nodes as visited
% record actions taken at each step (visited nodes at each step)
action_order = NaN(1, num_actions); 

Adj = adjacency(Ga); % adjacency matrix to figure out source nodes (current fringe)
Ga_unvisited = 1:num_actions; % list of unvisited nodes in Ga. run until empty.

curr_step = 1; 
while ~isempty(Ga_unvisited)
    % find frindge: those unvisited nodes with 0 column sum (no incoming edges)
    fringe = intersect(find(sum(Adj)==0), Ga_unvisited); 
    fringe_edges = []; 

    % create fringe graph with NTOC edges
    for g_ind = 1:numel(fringe)-1
        % for each action node, check if any og-nodes (nodes in 
        % original graph) in that action overlap with any other action node in
        % fringe. only need check other remaining nodes in fringe.
        for g_ind2 = g_ind+1:numel(fringe)
            % check if any nodes in action overlap 
            if ~isempty(intersect(Ga_n2a(fringe(g_ind), [2 3]), ...
                    Ga_n2a(fringe(g_ind2), [2 3])))
                % if there's a common node, add NTOC edge
                fringe_edges = [fringe_edges; g_ind g_ind2]; 
            end
        end
    end  

    % fringe graph with NTOC (undirected) edges [if no edges, then fine.]
    if isempty(fringe_edges)
        Gfr = graph([], [], 1, numel(fringe)); 
    else
        Gfr = graph(fringe_edges(:, 1), fringe_edges(:,2), 1, numel(fringe)); 
    end
    if plot_on 
        figure(); 
        plot(Gfr);
        title(sprintf('fringe graph%d', curr_step))
    end

    %% 3b. find max WIS (or approximate). We will take the best of WGMIN and WGMAX
    % ref: Sakai et al 2003 https://www.sciencedirect.com/science/article/pii/S0166218X02002056
    % calculate weights for each node using heuristic
    % heuristic 1 = 1 + # hops left for message that action is sending (if this
    % message path branches out in a subtree, sum all hops (edges) in the entire subtree)
    heur1 = @(node) 1+numel(nearest(Ga, fringe(node), Inf, 'Direction', 'outgoing')); % count all descendents
    fr_weights = arrayfun(heur1, 1:numel(fringe)); % fringe weights
    %% WGMIN greedy algorithm. 
    % input: graph with weighted nodes. output: greedily-selected maximal
    % independent set; using heuristic of weight(node)/(degree(node)+1)
    WGMIN_heur = @(node) fr_weights(node)./(degree(Gfr, node)+1); 
    WGMIN_weights = arrayfun(WGMIN_heur, 1:numel(fringe)); 
    % pick off valid nodes with highest weights until no more valid nodes
    WGMIN_set = []; 
    valid_nodes = 1:numel(fringe);
    while ~isempty(valid_nodes)
        [~, ni] = max(WGMIN_weights(valid_nodes));
        WGMIN_set = [WGMIN_set, valid_nodes(ni)]; % add chosen node to WGMIN_set
        nc = valid_nodes(ni); % chosen node (node # in Gfr) 
        % remove this chosen node, as well as all neighbors (Gfr with NTOC edges)
        % neighbors of selected node 
        neighs = neighbors(Gfr, nc)'; % row vector
        valid_nodes([ni, find(ismember(valid_nodes, neighs))]) = [];
    end

    %% WGMAX greedy algorithm
    % input: graph with weighted nodes. output: greedily-(remaining) maximal
    % independent set; using heuristic of
    % weight(node)/(deg(node)*(deg(node+1)) and removing nodes with smallest
    % heuristic value until no edges remain 

    WGMAX_heur = @(node) fr_weights(node)./(degree(Gfr, node)*(degree(Gfr,node)+1)); 
    WGMAX_weights = arrayfun(WGMAX_heur, 1:numel(fringe));

    WGMAX_set = 1:numel(fringe); % we'll take what node indices remains after deletions

    while ~isempty(subgraph(Gfr, WGMAX_set).Edges) % check if any remaining edges
        [~, ni] = min(WGMAX_weights(WGMAX_set)); % ni = index in WGMAX_set to be deleted
        % remove this chosen node
        WGMAX_set(ni) = [];
    end

    % take the bigger weighted (in Gfr weights, not greedy alg heuristic weights)
    % of WGMIN_set and WGMAX_set as the weighted independent set (WIS)
    if sum(fr_weights(WGMIN_set)) > sum(fr_weights(WGMAX_set))
        WIS = WGMIN_set; 
    else 
        WIS = WGMAX_set; 
    end
    %% 3c. mark WIS nodes as visited in action graph 
    Adj(fringe(WIS), :) = 0; % update adjacency matrix by clearing all marked nodes.
    Ga_unvisited(ismember(Ga_unvisited, fringe(WIS))) = []; 
    action_order(fringe(WIS)) = curr_step;
    
    curr_step = curr_step+1; 
end

%% 4. Scheduler
num_nodes0 = max(Ga_n2a(:, 2:3), [], 'all'); % oG nodes, involved in any message passing
node_schedules = cell(num_nodes0, 1); % each ith elt is node i's action schedule
msg_arrivals = nan(num_nodes0, num_nodes0); % each ith row = is node i's message's arrival time at nodes 1:num_nodes(G0)
% note if not all nodes have messages in m_paths, some msg_arrivals rows =
% nans, except the diagonal
msg_arrivals(logical(eye(size(msg_arrivals)))) = 0; % starting node already has message
action_schedule = zeros(num_actions, 5); 
act_ind = 0; % index of actions to make master schedule of actions
for stepi = 1:max(action_order)
    % batches of actions executed in each step 
    for act = find(action_order==stepi)
        act_ind = act_ind+1; 
        % for each action, check oG nodes involved 
        startnode = Ga_n2a(act, 2); 
        endnode = Ga_n2a(act, 3); 
        msg = Ga_n2a(act, 4); 
        
        % When adding a new action (send message i from j->k) action starts
        % at max(latest ending time of j, latest ending time of k) for both
        % nodes.
        nt = [node_schedules{startnode}; node_schedules{endnode}];
        if isempty(nt)
            tstart = 0; 
        else 
            tstart = max(nt(:, 2));
        end
        tstop = tstart + travel_time(startnode, endnode, msg); % time to execute action 
        node_schedules{startnode} = [node_schedules{startnode}; 
            tstart, tstop, msg];
        node_schedules{endnode} = [node_schedules{endnode}; 
            tstart, tstop, msg];
        action_schedule(act_ind, :) = [tstart, tstop, startnode, endnode, msg]; 
        msg_arrivals(msg, endnode) = tstop; 
    end
end
msg_arrivals;
% % MPT = max(msg_arrivals, [], 'all'); % message propagation time

