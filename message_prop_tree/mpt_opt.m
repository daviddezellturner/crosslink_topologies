%% brute force calculate the "optimal" root message passing time (MPT)
% Maryann Rui 01/14/2022 

% 1. choose a root node, form shortest path tree(based on distances or distances^2
% or other heuristic-- maybe a modified Dijkstra's with penalty for depth of tree?)
% 2. leaf nodes send message up to the root; then the conglomerate message is
% passed back down the tree (no need to send messages if a child node has
% already seen that message)

% 
% Notes:
% Each node can only receive one message at a time, and can only transmit
% one message at a time (to one other node). 
%% sample graph
%% generate graph and adjacency matrix
generate_new = false; 
if generate_new
    num_nodes = 8; % number of nodes
    A = max((randn(num_nodes) + 0.5),  0); % random undirected adjancency matrix, 
    % 0 = not connected, any positive value is distance btwn nodes. 
    % A(i,j) = d(i,j) [unfortunately there is a discontinuity where 0 =
    % infinitely far away but 0.1 is very very close]
    A = A - A.*eye(size(A)); 
    A = A + A'; 
end
A_d = A; % adjacency matrix based on d.
G_d = graph(A_d); % graph based on A
% A_d2 = A.^10; % adjacency matrix based on d^2. 
% G_d2 = graph(A_d2);

%% 
% RESULTS: d^2-shortest trees doesnt really change anything significantly.
shortest_tree_heuristic = 'd'; % 'd', 'd2'
% Create shortest path tree and calculate message passing times to the root
% and back down to the leaves
% -> see [nodes_tau, propagation_time_to_root] = calc_messageproptimes(Gd, root, m0)
num_nodes = size(G_d.Nodes, 1);
m0 = 1; % assume unit message sizes for now (each node has a size m0 message that needs to be transmitted to all other nodes)

% now brute force optimize for propagation_time_to_root and
% propagation_time_from_root
root_latencies = zeros(1, num_nodes); 
total_latencies = zeros(1, num_nodes); 
for nodei = 1:num_nodes
    [propagation_time_to_root, propagation_time_from_root, ~, ~] = calc_messageproptimes(G_d, nodei, m0, shortest_tree_heuristic); 
    root_latencies(nodei) = propagation_time_to_root; 
    total_latencies(nodei) = propagation_time_from_root + propagation_time_to_root; 
end
root_latencies
total_latencies
[~, t_opt_node] = min(root_latencies)
[~, t_tot_opt_node] = min(total_latencies)
% so the optimal root node is actually 4, which does not coincide with any
% centrality based measure that we've tried so far. 
[~, ~, nodes_tau_up] = calc_messageproptimes(G_d, t_opt_node , m0, shortest_tree_heuristic);

    
    
%% plot results: tree overlaid on original graph with latencies
figure(); 
H0 = plot(G_d);   % original plot
% tree plot
Gt = shortestpathtree(G_d, t_tot_opt_node);
highlight(H0, Gt); % show root's shortest path tree on orig graph
highlight(H0, t_opt_node, 'NodeColor', '#EDB120'); % highlight the root yellow
labeledge(H0, Gt.Edges.EndNodes(:, 1), Gt.Edges.EndNodes(:, 2),...
    round(Gt.Edges.Weight, 3)); % label edges involved in sptree
labelnode(H0, 1:num_nodes, round(nodes_tau_up,5));
H0.NodeLabelColor='#D95319';
H0.NodeFontSize=11;
title('Message propagation tree with latencies');% , d2-trees')
% saveas(gcf, 'sample_message_prop_tree.png');
%% plot results: original graph with (up) MPT (or latency) for each node as root
figure();
H1 = plot(G_d); 
labelnode(H1, 1:num_nodes, round(root_latencies, 3));
title('Upward propagation latencies for each node if it were the root'); %, using d2-trees')
highlight(H1, t_tot_opt_node, 'NodeColor', '#EDB120');
H1.NodeFontSize=11;
% saveas(gcf, 'sample_upwd_prop_times.png');
%% plot results: original graph with total MPT (up + down) for each node as root 
figure();
H2 = plot(G_d); 
labelnode(H2, 1:num_nodes, round(total_latencies, 3));
title('Total message propagation times for each node if it were the root'); %, using d2-trees')
highlight(H2, t_tot_opt_node, 'NodeColor', '#EDB120');
H2.NodeFontSize=11;
% saveas(gcf, 'sample_total_prop_times.png');
%% plot: graph with node labels; 
figure(); 
H3 = plot(G_d); 
H3.NodeFontSize=11;
title('Graph with node labels'); 
% saveas(gcf, 'sample_graph_labeled.png');
% MPT = message propagation time
