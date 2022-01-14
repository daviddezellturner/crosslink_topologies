%% INPUT: a directed graph, a source node A, and a target node B.
%% OUTPUT: the shortest path from A to B, and the associated cost.
%% if no such path exists, best_path will be [] and cost will be -1. 

% the pathfinder function assumes that each node already has a complete 
% picture of the network. beginning with such a picture, it takes in a 
% node-pair (A,B) and determines the shortest path from A to B at a given 
% timeslot, given a set of weights. (this is pretty much just djikstra's
% algorithm adapted to our input format! the main advantage of having it
% implemented is to more easily track our metrics dynamically if we want..)

% ultimately, we will feed this function into our genetic algorithm (after
% routing) to optimize the algorithm over many timesteps.

% IMPORTANT NOTE: we assume that a zero weight in the input graph implies 
% nonexistence of edge from row node to column node. also, input graph is
% directed!

% A and B are the indices of source and target nodes.
function [best_path, cost] = pathfinder(graph, A, B)

% mark all nodes unvisited. create a of all unvisited nodes.

num_nodes = size(graph,1);
current_node = A;
unvisited_set = zeros(1,num_nodes);
for node = 1:num_nodes
    unvisited_set(node)=node;
end

% initialize distances. (begin with zero for node A, and infinity for
% everything else. we will represent infinity via -1 since this is not a
% possible cost value (we are constraining ourselves to positive values, at
% least for now). also initialize paths.

costs = zeros(1,num_nodes);
paths = cell(1,num_nodes);
for node = 1:num_nodes
    paths(1,node) = {[current_node]};
    if node == current_node
        costs(node)=0;
    else
        costs(node)=-1;
    end
end

% if destination node has been marked visited, or if smallest tentative
% distance in unvisited set is infinity, then stop. otherwise, select 
% unvisited node marked w/ smallest tentative distance, set as new current 
% node, and continue to path-find.

reached_destination = 0;
unreachable = 0;

cost = -1;
best_path = [];

while ~reached_destination && ~unreachable

    % for current node, consider all unvisited neighbors, and calculate
    % tentative distances through current node. compare newly computed distance
    % to current, and assign smaller one. then mark current node as visited.
    
    for node = 1:num_nodes
        if node ~= current_node
            if graph(current_node, node) ~= 0
                alternate_cost = costs(current_node) + graph(current_node, node);
                if costs(node) == -1
                    costs(node) = alternate_cost;
                    paths(node) = {[paths{current_node} node]};
                else
                    if min([costs(node), alternate_cost]) == alternate_cost
                        costs(node) = alternate_cost;
                        paths(node) = {[paths{current_node} node]};
                    end
                end
            end
        end
    end
    
    unvisited_set(unvisited_set==current_node) = [];

    % check if we've finished path-finding
    if ~ismember(B, unvisited_set)
        cost = costs(B);
        best_path = paths(B);
        reached_destination = 1;
    end
    
    % find lowest-cost neighbor from which to continue search
    smallest_unvisited_cost = -1;
    smallest_unvisited_node = current_node;
    for unvisited_node = unvisited_set
        if smallest_unvisited_cost == -1
            smallest_unvisited_cost = costs(unvisited_node);
            smallest_unvisited_node = unvisited_node;
        elseif costs(unvisited_node) ~= -1
            if min(smallest_unvisited_cost, costs(unvisited_node)) == costs(unvisited_node)
                smallest_unvisited_cost = costs(unvisited_node);
                smallest_unvisited_node = unvisited_node;
            end
        end
    end
    current_node = smallest_unvisited_node;

    % check for unreachability
    if smallest_unvisited_cost == -1
        unreachable = 1;
    end

end