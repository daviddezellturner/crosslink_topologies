%% Simulates a link-state routing protocol
% each node will gain a picture of the whole network
% (cost of pings is probably negligible compared to message-sending, so this
% isn't super important, but might be nice to have/think about later)

function network = linkstate(graph)

% In our simulation, each node starts off only having information about its
% direct neighbors (e.g. by pinging). Nodes build up their map via LSAs.

nodes = [];
num_nodes = size(graph,1);
sequence_nums = ones(1,num_nodes);
known_info = cell(1,num_nodes);

% Contents of a link state announcement (LSA): origin adress, sequence
% number, (neighbor1, cost1), (neighbor2, cost2), ...., (neighborn, costn)

for node = 1:num_nodes
    adjacencies = graph(node,:);
    known_info(node) = [adjacencies sequence_nums(node)];
end

% process might finish sooner in practice -- optimze
for timestep = 1:num_nodes
    for node = 1:num_nodes
        for node2 = 1:num_nodes
            if graph(node,node2) ~= 0
                known_info(node2) = [known_info{node2} known_info{node}];
            end
        end
    end
end

% re-broadcast with upated sequence number every time a new neighbor/link is added or deleted
% in practice, rather than doing any of the above, maybe just run the whole
% procedure every x seconds

% placeholder (result of constructing from adjacencies will match)
network = graph;