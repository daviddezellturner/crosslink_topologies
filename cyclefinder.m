% To test: how will traversing a greedy cycle (2x?) compare with a
% bottom-up tree approach? In terms of local vs. global optimization?
% can we trade off b/t greedy optimization and some amount of
% "look-forward"? dynamic programming, back-tracking, implementation of
% existing algos, ...

%% INPUT: a directed graph and a "root node."
%% OUTPUT: a greedy cycle. if no cycle exists, return path [] and a -1 cost
% may also return no path if greediness causes failure - think abt how to fix

% A and B are the indices of source and target nodes.
function [greedy_cycle, cost] = cyclefinder(graph, root)

    num_nodes = size(graph,1);
    current_node = root;
    unvisited_set = zeros(1,num_nodes);
    for node = 1:num_nodes
        unvisited_set(node)=node;
    end
    greedy_cycle = [current_node];

    unvisited_set(unvisited_set==current_node) = [];
    
    closest_neighbor = current_node;
    cost = 0;
    while length(greedy_cycle) < num_nodes
        min_dist = -1;
        for node = unvisited_set
            if graph(current_node,node) ~= 0
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
        if min_dist == -1
            greedy_cycle = [];
            cost = -1;
            return
        else
            greedy_cycle = [greedy_cycle closest_neighbor];
            unvisited_set(unvisited_set==closest_neighbor) = [];
            cost = min_dist + cost;
            current_node = closest_neighbor;
        end
    end
    if graph(current_node, root) ~= 0
        greedy_cycle = [greedy_cycle root];
        cost = cost + graph(current_node, root);
    else
        greedy_cycle = [];
        cost = -1;
        % should change this to try again by starting w/ second-closest
        % neighbor!!!
    end
end