%% INPUT: a directed graph and a "root node."
%% OUTPUT: a greedy cycle (hamiltonian or near-hamiltonian). 
%% if no cycle exists, return path [] and a -1 cost.

% To test: how will traversing a greedy cycle compare with a
% bottom-up tree approach, with channel-sharing?

% added in a heuristic function as in a* ("distance to target" is probably
% not the right kind of heuristic function to use in this case, b/c we are
% trying to create a cycle which hits all nodes -- something along the
% lines of "value having lots of neighbors" could be better)

% allow returning non-hamiltonian cycles (i.e. some revisiting
% allowed) if this ends up being necessary and/or "better" by some metric

% (might also want to adopt idea of heuristics to the tree methodology)

function [greedy_cycle, greedy_cost] = cyclefinder_backtrack_heuristic(graph, root)

    num_nodes = size(graph,1);
    greedy_cycle=[root];
    greedy_cost=-1;
    found_cycle=false;

    visited = zeros(1,num_nodes);
    revisits = zeros(1,num_nodes);
    revisit_cap = 2;

    function find_cycle(node)

        %greedy_cycle
        %revisits
        
%         if length(greedy_cycle) > num_nodes*3/2
% 
%             visited(greedy_cycle(length(greedy_cycle))) = 1;
%             greedy_cycle(length(greedy_cycle)) = [];
%             visited(greedy_cycle(length(greedy_cycle))) = 1;
%             greedy_cycle(length(greedy_cycle)) = [];
% 
%         end

        % check if we've found a cycle
        if node == root & length(greedy_cycle) >= num_nodes+1
            missing=false;
            for n=1:num_nodes
                if ~ismember(n,greedy_cycle)
                    missing=true;
                    break
                end
            end
            if ~missing
                found_cycle=true;
                return
            end
        end

        % otherwise check if we can close a cycle
        if node ~= root & length(greedy_cycle) >= num_nodes
            missing=false;
            for n=1:num_nodes
                if ~ismember(n,greedy_cycle)
                    missing=true;
                    break
                end
            end
            if ~missing
                if graph(greedy_cycle(length(greedy_cycle)),root) ~= 0
                    greedy_cycle = [greedy_cycle root];
                    found_cycle = true;
                    return
                end
            end
        end

        % else, continue to explore

        neighbors_min_to_max = [];
        for i = 1:num_nodes
            current_min_cost = -1;
            current_min_index = -1;
            avg = mean(graph);
            avg = avg(1);
            for j = 1:num_nodes
                num_neighbors = length(graph(j,:));
                cost_function = 0.5*graph(node,j)/avg+0.5*num_neighbors/num_nodes;
                if ~ismember(j,neighbors_min_to_max) & graph(node,j) ~=0 ...
                        & (current_min_cost == -1 | cost_function <= current_min_cost)
                    current_min_cost = cost_function;
                    current_min_index = j;
                end
            end
            if current_min_cost ~= -1
                neighbors_min_to_max = [neighbors_min_to_max current_min_index];
            end
        end

        % explore via recursion on neighbors
        for neighbor = neighbors_min_to_max

            
            if visited(neighbor) == 0 & revisits(neighbor) < revisit_cap
                visited(neighbor)=1;
                revisits(neighbor)=revisits(neighbor)+1;
                greedy_cycle = [greedy_cycle neighbor];
                find_cycle(neighbor);
                if found_cycle
                    return
                end
                visited(neighbor)=0;
                revisits(neighbor)=revisits(neighbor)-1;
                greedy_cycle(length(greedy_cycle)) = [];
            end
        end

        % allowing revisits
        if length(greedy_cycle) >= num_nodes*4/5

            visited=zeros(1,num_nodes);
            visited(node)=1;

            for neighbor = neighbors_min_to_max
                if visited(neighbor) == 0 & revisits(neighbor) < revisit_cap
                    visited(neighbor)=1;
                    revisits(neighbor)=revisits(neighbor)+1;
                    greedy_cycle = [greedy_cycle neighbor];
                    find_cycle(neighbor);
                    if found_cycle
                        return
                    end
                    visited(neighbor)=0;
                    revisits(neighbor)=revisits(neighbor)-1;
                    greedy_cycle(length(greedy_cycle)) = [];
                end
            end
            
        end

    end

    % recursive cycle computation
    find_cycle(root);

    % final cost calculation
    % this isn't accurate anymore -- see computation in testing
    if found_cycle
        greedy_cost=0;
        for node = 1:length(greedy_cycle)-1
            greedy_cost=greedy_cost+graph(greedy_cycle(node),greedy_cycle(node+1));
        end
    end

end