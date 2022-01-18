%% INPUT: a directed graph and a "root node."
%% OUTPUT: a greedy cycle (hamiltonian or near-hamiltonian). 
%% if no cycle exists, return path [] and a -1 cost.

% To test: how will traversing a greedy cycle compare with a
% bottom-up tree approach? In terms of local vs. global optimization?
% can we trade off b/t greedy optimization and some amount of
% "look-forward"? dynamic programming, back-tracking, implementation of
% existing algos, ...

% for now this is just a greedy bfs; could try to add in some kind of
% heuristic function to be more like a* ("distance to target" is probably
% Not the right kind of heuristic function to use in this case, b/c we are
% trying to create a cycle which hits All nodes -- something along the
% lines of "value having lots of neighbors" could be better)

% could also explore returning non-hamiltonian cycles (i.e. some revisiting
% allowed) if this ends up being "better" by some metric

function [greedy_cycle, cost] = cyclefinder_backtrack(graph, root)

    num_nodes = size(graph,1);
    current_node = root;
    unvisited_set = zeros(1,num_nodes);
    for node = 1:num_nodes
        unvisited_set(node)=node;
    end

    greedy_cycle = [current_node];

    unvisited_set(unvisited_set==current_node) = [];

    closest_neighbor = current_node;
    avoiding = [];
    cost = 0;
    while length(greedy_cycle) < num_nodes
        greedy_cycle
        min_dist = -1;
        for node = unvisited_set
            if graph(current_node,node) ~= 0 & isempty(avoiding(avoiding==current_node))
                if min_dist == -1
                    min_dist = graph(current_node,node);
                    closest_neighbor = node;
                else
                    if min(min_dist, graph(current_node,node)) == graph(current_node,node)
                        min_dist = graph(current_node,node);
                        closest_neighbor = node;
                    end
                end
            end
        end

        % BACKTRACK
        if min_dist == -1
            avoiding = [avoiding current_node];
            avoiding
            prev_layer = length(greedy_cycle);
            if prev_layer == 1
                cost = 0;
            else
                cost = cost - graph(greedy_cycle(prev_layer),current_node);
            end
            greedy_cycle(greedy_cycle==current_node) = [];
            % backtrack another layer
            if length(avoiding) + length(greedy_cycle) == num_nodes
                if length(greedy_cycle) == 0
                    cost = -1;
                    return
                else
                    prev_node = greedy_cycle(length(greedy_cycle));
                    greedy_cycle(greedy_cycle==current_node) = [];
                end
            end
        else
            greedy_cycle = [greedy_cycle closest_neighbor];
            unvisited_set(unvisited_set==closest_neighbor) = [];
            cost = min_dist + cost;
            current_node = closest_neighbor;
        end
    end


    % final path (TODO add in special backtracking step for case where we
    % can't close the loop)
    if graph(current_node, root) ~= 0
        greedy_cycle = [greedy_cycle root];
        cost = cost + graph(current_node, root);
    else
        greedy_cycle = [];
        cost = -1;
    end
end
