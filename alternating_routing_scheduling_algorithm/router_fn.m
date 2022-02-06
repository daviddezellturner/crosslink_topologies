function [best, final, alg_MPTs, alg_deltas] = router_fn(A0, source_nodes, dest_nodes, travel_time, SIM_ANNEAL, verbose)
% Maryann Rui 
% January 2022 

% RETURNS: best and final routing schedules of algorithm
%   - alg_MPTs, alg_deltas = mpts and normalized change in edge-MPT
%   estimates per iteration of the algorithm

%% 1. Router
% Given: 
% > Graph with estimated edge-mpt's in a 3d matrix with indices (startnode, endnode, message)
% > Source nodes
% > Destination nodes (same set for all source nodes)
%% 1a. initialize graph and estimated edge-mpt times
num_nodes_G = size(A0, 1); 
G0 = graph(A0);

% Initial estimates = simply edge travel time for given message 
[x,y,z] = ndgrid(1:num_nodes_G, 1:num_nodes_G, 1:num_nodes_G);

est_mpt = arrayfun(travel_time, x,y,z);
est_mpt(est_mpt==Inf) = 0;
%% 2b. Find message routes via shortest path tree w/ edge-mpt estimates

MAX_NUM_ITERS = 200; 
iter = 0; 
% keep track of best MPT so far (since we dont necessarily converge to min)
best.MPT = Inf; 

alg_MPTs = zeros(1, MAX_NUM_ITERS); % record MPT for each iteration
alg_deltas = zeros(1, MAX_NUM_ITERS); % record est_change [normalized by 
% previous edge estimate matrix norm(2-norm)] (in edge estimates) for each iteration
est_change = Inf; 

while est_change > 0.0001 && (iter < MAX_NUM_ITERS) 
    iter = iter + 1; 
    tmp_paths = cell(1, numel(source_nodes)); 
    for ind = 1:numel(source_nodes)
        ns = source_nodes(ind);
        G_est_i = digraph(est_mpt(:, :, ns)); % estimated edge times on graph, for message i
        Gi = shortestpathtree(G_est_i, ns, dest_nodes);
        tmp_paths{ind} = Gi; 
        if PLOT_ITERS && ind==1
            highlight(H1, Gi); 
        end
    end
    m_paths = containers.Map(num2cell(source_nodes), tmp_paths);

    %% 2. Scheduler-Router interface
    plot_scheduler = false; 
    [node_schedules, msg_arrivals, Ga_n2a, action_schedule] = scheduler_fn(m_paths, travel_time, plot_scheduler);
    MPT = max(msg_arrivals, [], 'all'); % message propagation time
    if MPT < best.MPT
        best.MPT = MPT; 
        best.node_schedules = node_schedules; 
        best.msg_arrivals = msg_arrivals;  
        best.m_paths = m_paths;  
        best.action_schedule = action_schedule;
        best.iter = iter; 
    end
    if verbose
        fprintf('Iter %d: MPT = %0.2f s\n ', iter, MPT); 
    end

    %% 3. Updating estimated routing graph edge times
    % edge-MPT = (time it takes for a message i to cross j->k after arriving at
    % j, taking into account waiting in the queue at j plus the travel time
    % j->k)
    edge_mpt = @(j, k, i) msg_arrivals(i, k) - msg_arrivals(i,j); % message i, j->k
    est_mpt_prev = est_mpt; 
    
    for ind = 1:size(Ga_n2a, 1)
        action = Ga_n2a(ind, 2:4); % [startnode, endnode, msg]
        % update edge weight
        est_mpt(action(1),action(2),action(3)) = edge_mpt(action(1), action(2), action(3));
    end
    % est_mpt is updated
    est_change = norm(est_mpt_prev(:) - est_mpt(:))/norm(est_mpt_prev(:));
    if verbose
        fprintf('Normalized change in edge-mpt estimates: %.4g \n', est_change);
    end

    if PLOT_ITERS
        figure(); 
        H = plot(G0); 
        msg_to_plot = 1; 
        label_msg_arrivals(H, msg_arrivals, msg_to_plot);
        title(sprintf('Iter %d: Graph with Message %d Arrival Times', iter, msg_to_plot))
    end
    
    alg_MPTs(iter) = MPT; 
    alg_deltas(iter) = est_change; 
    
    %% 4. Simulated Annealing- Add noise to est_mpt (estimated edge-message propagation times) 
    if SIM_ANNEAL 
        % only add noise to nonzero elements to ensure 0 => no connection
        est_mpt = est_mpt + randn(size(est_mpt))*(sqrt(max(est_mpt(:)))/6).*logical(est_mpt);
        % clip negative values to 0.1
        est_mpt(est_mpt<0) = 0.1;
    end
    
end
% save details of final iteration
final.MPT = MPT; 
final.node_schedules = node_schedules; 
final.msg_arrivals = msg_arrivals;  
final.m_paths = m_paths;  
final.action_schedule = action_schedule; 
final.iter = iter; 
if verbose
    fprintf('The best MPT occured at iter %d, with MPT = %.4g s\n', best.iter, best.MPT); 
end
iters = 1:iter; 
alg_MPTs = alg_MPTs(iters); 
alg_deltas = alg_deltas(iters); 

end