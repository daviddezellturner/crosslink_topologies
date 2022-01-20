%% INPUT: a directed graph and a "root node."
%% OUTPUT: a greedy cycle (hamiltonian or near-hamiltonian). 
%% if no cycle exists, return path [] and a -1 cost.

% To test: how will traversing a greedy cycle compare with a
% bottom-up tree approach? In terms of local vs. global optimization?

% for now this is just a greedy bfs; could try to add in some kind of
% heuristic function to be more like a* ("distance to target" is probably
% not the right kind of heuristic function to use in this case, b/c we are
% trying to create a cycle which hits all nodes -- something along the
% lines of "value having lots of neighbors" could be better)

% could also explore returning non-hamiltonian cycles (i.e. some revisiting
% allowed) if this ends up being "better" by some metric

% (might also want to adopt idea of heuristics to the tree methodology)

function [greedy_cycle, greedy_cost] = cyclefinder_backtrack(graph, root)

    num_nodes = size(graph,1);
    greedy_cycle=[root];
    greedy_cost=-1;
    found_cycle=false;

    visited = zeros(1,num_nodes);

    function find_cycle(node)

        % check if we've found a cycle
        if node == root & length(greedy_cycle) == num_nodes+1
            found_cycle=true;
            return            
        end

        % if we haven't found a cycle, continue to explore

        % computing ordered list of neighbors to prioritize
        % TODO: incorporate heuristic function here!
        neighbors_min_to_max = [];
        for i = 1:num_nodes
            current_min_cost = -1;
            current_min_index = -1;
            for j = 1:num_nodes
                if ~ismember(j,neighbors_min_to_max) & graph(node,j) ~=0 ...
                        & (current_min_cost == -1 | graph(node,j) <= current_min_cost)
                    current_min_cost = graph(node,j);
                    current_min_index = j;
                end
            end
            if current_min_cost ~= -1
                neighbors_min_to_max = [neighbors_min_to_max current_min_index];
            end
        end

        % explore via recursion on neighbors
        for neighbor = neighbors_min_to_max
            if visited(neighbor) == 0
                visited(neighbor)=1;
                greedy_cycle = [greedy_cycle neighbor];
                find_cycle(neighbor);
                if found_cycle
                    return
                end
                visited(neighbor)=0;
                greedy_cycle(length(greedy_cycle)) = [];
            end
        end

    end

    % recursive cycle computation
    find_cycle(root);

    % final cost calculation
    if found_cycle
        greedy_cost=0;
        for node = 1:num_nodes-1
            greedy_cost=greedy_cost+graph(greedy_cycle(node),greedy_cycle(node+1));
        end
        greedy_cost=greedy_cost+graph(greedy_cycle(num_nodes),greedy_cycle(1));
    end

end
