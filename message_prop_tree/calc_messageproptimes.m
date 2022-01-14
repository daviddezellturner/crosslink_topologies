function [propagation_time_to_root, propagation_time_from_root, nodes_tau_up, nodes_tau_down] = calc_messageproptimes(Gd, root, m0, sptree_heuristic)
% %%%% INPUTS %%%%%%
% Gd = graph with distances as edge weights 
% root = a node in the graph to be the root of the tree
% m0 = size of each node's message, assumed uniform for now. (in the future
% can make it a node property)
% sptree_heuristic = 'd' or 'd2' for shortest tree criteria. tree for message passing.
% %%%%  RETURNS %%%%%
% nodes_tau = propagation time of a message leaves up to given node; 
% we're interested in propagation_time_to_root = nodes_tau(root);

% use 'd' or 'd2' shortest path tree as heuristic for best tree from root
% for message passing. 
if nargin==4
    if isequal(sptree_heuristic, 'd2')
        Ad2 = adjacency(Gd, 'weighted').^2; 
        Gd2 = graph(Ad2);
        Gt = shortestpathtree(Gd2, root);
    else % use 'd' shortest tree
        Gt = shortestpathtree(Gd, root);
    end
else 
    Gt = shortestpathtree(Gd, root); % d-shortest tree is default
end

% Calculate each node's level in the tree: start from root (level 0)
node_levels = calc_levels(Gt, root);  % 1xn
% Assign weights to nodes based on size of all their children's messages
node_weights = calc_node_weights(Gt, node_levels, m0); % nx1

% Calculate propagation time of a message from leaves to root.
nodes_tau_up = calc_upwd_propagation_times(Gt, Gd, node_weights, node_levels) 
propagation_time_to_root = nodes_tau_up(root);

nodes_tau_down = calc_dwdprop_times(Gt, Gd, node_weights, node_levels)
propagation_time_from_root = nodes_tau_down(root);
%% helper functions

function levels = calc_levels(Gt, root)
    % calculates which level each node is on. root is on level 0.
    % Gt should be a tree. 
    % returns row vector levels (1xnum_nodes) of node levels
    level=0; 
    num_nodes = size(Gt.Nodes, 1);
    levels = zeros(1, num_nodes); 
    level_nodes = root; 
    while ~isempty(level_nodes)
        levels(level_nodes) = level; 
        level = level + 1; 
        new_level = []; 
        for lp = level_nodes' % level parent
            % make children the new level_nodes
            new_level = [new_level; successors(Gt, lp)];
        end
        level_nodes = new_level;
    end
end

    
function node_weights = calc_node_weights(Gt, node_levels, m0)
    % WEIGHTS: moving from bottom level up, assign nodes weights based on 
    % how many children they have, +m0 for themselves
    % input: Gt = directed tree, node_levels = 1xnum_nodes vector of which
    % tree level nodes are on (root = level 0)
    % m0 = size of each node's message.
    num_nodes = size(Gt.Nodes,1); 
    node_weights = zeros(num_nodes,1); 
    num_levels = max(node_levels);
    for levelnum = num_levels:-1:0
        nodes = find(node_levels==levelnum); % nodes at level
        % for each node, +1 to itself, +1 to parent
        node_weights(nodes) = node_weights(nodes)+m0; % self weight
        for ni = nodes % make sure it's row vector
            pni = predecessors(Gt, ni);
            node_weights(pni) = node_weights(pni)+node_weights(ni);
        end
    end
end

function [nodes_tau, edge_times] = calc_upwd_propagation_times(Gt, Gd, node_weights, node_levels)
    % inputs: Gt = directed tree for message prop [disregard edge values]
    % Gd = original graph (needed for distances)
    % node_weights = sum of node's weight + all childrens' original message
    % sizes
    % node_levels (1xnum_nodes) level of each node in the tree (root = 0)
    % interested in nodes_tau(root node) 
    
    num_nodes = size(Gt.Nodes,1); 
    Adj = adjacency(Gd, 'weighted'); % adjacency matrix with distances
    % Calculate time to transmit signal of size node_weights(i) from node i
    % to node j, given distance Adj(i,j) between nodes i, j. 
    edge_times = Adj.^2 .* repmat(node_weights, 1, num_nodes);
    
    num_levels = max(node_levels);
    nodes_tau = zeros(num_nodes,1); % in particular, leaf tau's are 0
    for levelk = (num_levels-1):-1:0
         nodes = find(node_levels==levelk); % nodes at level
         for nodei = nodes
            % tau_i = max_{k = child(i)} (tau_k + t_{ki}) 
            childi = successors(Gt, nodei);
            if ~isempty(childi) % if not a leaf node
                % order children based on tau's
                childtaus = [nodes_tau(childi), childi];
                childtaus = sortrows(childtaus); % sort based on first column
                zi = zeros(1, numel(childi)); 
                zi(1) = childtaus(1,1) + edge_times(childtaus(1,2), nodei);
                for indi = 2:numel(childi)
                    zi(indi) = max(zi(indi-1), childtaus(indi, 1)) + edge_times(childtaus(indi, 2), nodei); 
                end
                nodes_tau(nodei) = zi(numel(childi));
            end
         end
    end
end

function [nodes_tau, edge_times] = calc_dwdprop_times(Gt, Gd, node_weights, node_levels)
    % node_weights nx1 
    
    % Here, nodes_tau(i) = time it takes for node i to pass the part of total message
    % that its children doesn't have, and so on, all the way down to the leaves
    % of its subtree. We are interested in nodes_tau(0)
    
    
    num_levels = max(node_levels); 
    % adjacency matrix with original distances from G_d
    Adj = adjacency(Gd, 'weighted'); 
    num_nodes = size(Gd.Nodes, 1); 
    % edge times are calculated based on a message of size (m0*num_nodes -
    % node_weight(j)) for sending a message from parent i to node j. 
    % it is directional.
    M = max(node_weights); % weight of all messages of the tree (should be m0*num_nodes) = node_weights(root)
    edge_times = Adj.^2 .* repmat((M - node_weights)', num_nodes, 1);

    % parent node can only broadcast to one child at a time
    
    % go from leaves up to root (level by level) and calculate tau's to
    % assign child-broadcasting order and calculate parent level's tau. 
    nodes_tau = zeros(num_nodes,1); % in particular, leaf tau's are 0
    for levelk = (num_levels-1):-1:0
         nodes = find(node_levels==levelk); % nodes at level
         for nodei = nodes
            childi = successors(Gt, nodei);
            if ~isempty(childi) % if not a leaf node
                % order children based on tau's in DESCENDING order
                childtaus = [nodes_tau(childi), childi]; %(tau, nodeid)
                childtaus = sortrows(childtaus, 'descend'); % sort based on first column
                zi = zeros(1, numel(childi)); % times it takes for children to finish sending parent's message down to its leaves
                runningt = 0; 
                for indi = 1:numel(childi)
                    runningt = runningt+edge_times(nodei, childtaus(indi, 2)); 
                    zi(indi) = runningt + childtaus(indi, 1); 
                end
                nodes_tau(nodei) = max(zi); 
            end
         end
    end
end

end
