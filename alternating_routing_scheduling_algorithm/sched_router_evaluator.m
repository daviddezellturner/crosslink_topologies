%% sched_router_evaluator.m 
% Maryann Rui
% January 2022 

% Runs router_fn.m over multiple time points of ephemerides 
% also plots output of router_fn (MPT-best and final routing schedules
% output by the algorithm) 
%% 
for test = "MMS" %["LunaNet", "ATrainCTrain", "MMS"]
%test = 'GPSWalker'
test = char(test);
switch test
    case 'LunaNet'
        S = load('./LunaNetGraphList.mat'); 
        GraphList=S.LunaNetGraphList; 
        timekeys = keys(GraphList);
    case 'ATrainCTrain'
        S = load('./ATrainCTrainGraphList.mat'); 
        GraphList = S.ATrainCTrainGraphList; 
        timekeys = keys(GraphList);
    case 'GPSWalker'
        S=load('./GPSWalkerApproxGraphList.mat');
        GraphList= S.GPSWalkerApprox;
        timekeys = keys(GraphList);
    case 'MMS'
        S=load('./MMSGraphList.mat');
        GraphList= S.MMSGraphList;
        timekeys = keys(GraphList);
end

% algorithm options
USE_DATA_RATE_TABLE = false; 
SIM_ANNEAL=false;

% post-processing options
PLOT_BASE= false; % plot basics: initial graph, alg output graph, best MPT graph, alg run info
SAVE_BASE = false; 

append='_x40';

SAVE_OUTPUTS=false;


%% run router-scheduler algorithm at multiple time points
numt = 50; 
tkeyinds = 1:50:numel(timekeys); 
final_MPT_overT = inf(1, numel(tkeyinds)); 
best_MPT_overT = inf(1, numel(tkeyinds)); 
T_vals = zeros(1, numel(tkeyinds)); 
for ind = 1:numel(tkeyinds)

    tkeyind = tkeyinds(ind);

    % get the right graph at given time
    T = timekeys{tkeyind}; 
    T_vals(ind) = T; 
    fprintf('T: %.3g s \n', T); 
    A0 = GraphList(T); % adj matrix at t = tkey
    num_nodes_G = size(A0, 1); 
    m_sizes = ones(1, num_nodes_G); 
    source_nodes = 1:num_nodes_G; 
    dest_nodes = 1:num_nodes_G;

    test_str = [test '_T' num2str(T) '_2'];

    G0 = graph(A0);

    if ~USE_DATA_RATE_TABLE
        % option 1: use exponential fit for data rate & travel_times
        a = 100.08; b = 0.667; 
        travel_time = @(j, k, i) m_sizes(i).*A0(j,k).^b./a; % travel time of message size m_i [m_sizes incorporated] along edge j->k
        dataratesource='rate_expfit';
    else
        % option 2: use slant range-data rate table for travel_times [see function below]
        travel_time = @(j,k,i) travel_time_table(j,k,i,A0, m_sizes);
        dataratesource='rate_table';
    end
    
    %% function - output interface
    verbose = false; 
    [best, final, alg_MPTs, alg_deltas] = router_fn(A0, source_nodes, dest_nodes, travel_time, SIM_ANNEAL, verbose); 
    fprintf('The best MPT occured at iter %d, with MPT = %.4g s\n', best.iter, best.MPT); 
    fprintf('The final MPT occured at iter %d, with MPT = %.4g s\n', final.iter, final.MPT); 
    final_MPT_overT(ind) = final.MPT; 
    best_MPT_overT(ind) = best.MPT; 
end
%% save data 
if SAVE_OUTPUTS
    save([test '_outputs' append '.mat'], 'best_MPT_overT','final_MPT_overT', 'T_vals'); 
    % % load('temp_save.mat');
    %%
    figure(); 
    plot(T_vals, best_MPT_overT); 
    xlabel('T_0 [s]'); ylabel(['Best MPT [s]']); 
    title([test ' Message Propagation Times over Ephemeris']); 
    saveas(gcf, [test '_mpt_over_time' append '.png']); 

end
end

%% 1. Plot algorithm run
if PLOT_BASE
    f0 = figure(); 
    H = plot(G0, 'EdgeLabel', round(G0.Edges.Weight,2)); 
    title([test ' Graph G_0 with slant ranges [km]']); 
    if SAVE_BASE
        saveas(f0, [test_str '_graph_dist.png']);
    end
    a = 100.08; b = 0.667; 
    edgetimesplt = A0.^b./a; % travel time of message size m_i [m_sizes incorporated] along edge j->k
    Gtime = graph(edgetimesplt); 
    f00 = figure(); 
    H = plot(Gtime, 'EdgeLabel', round(Gtime.Edges.Weight,3)); 
    title([test ' Graph G_0 with 1kb edge travel times [s]']); 
    if SAVE_BASE
        saveas(f00, [test_str '_graph_times.png']);
    end
    
    f1 = figure(); 
    plot(1:final.iter, alg_MPTs); 
    if SIM_ANNEAL
        title('Total Message Propagation Time (MPT) versus iteration, with Sim. Anneal'); 
    else
        title('Total Message Propagation Time (MPT) versus iteration'); 
    end
    xlabel('Iteration number'); 
    ylabel('MPT [s]'); 
    
    f2 = figure(); 
    plot(1:final.iter, alg_deltas);
    if SIM_ANNEAL
        title('Normalized change in edge-MPT estimate per iteration, with Sim. Anneal');
    else
        title('Normalized change in edge-MPT estimate per iteration');
    end
 
    xlabel('Iteration number'); 
    ylabel('Normalized change in edge-message propagation time'); 
    

    if SAVE_BASE 
        if SIM_ANNEAL
            saveas(f1, [test_str '_alg_run_mpt_sim_anneal_' dataratesource '.png']);
            saveas(f2, [test_str '_alg_run_deltas_sim_anneal_' dataratesource '.png']); 
        else 
            saveas(f1, [test_str '_alg_run_mpt_' dataratesource '.png']);
            saveas(f2, [test_str '_alg_run_deltas_' dataratesource '.png']); 
        end
    end
end

%% 2. Plot sample message path
if PLOT_BASE
    msg_to_plot = 2; 
    plot_things = {final, best}; 
    plot_strings = {'Final', 'Best'}; 
    % plot both: algorithm final output, best MPT seen in alg run

    for ind = 1:numel(plot_things)
        thing_to_plot = plot_things{ind}; 
        string_to_plot = plot_strings{ind};
        figure();
        H = plot_msg_path(G0, thing_to_plot.m_paths, msg_to_plot, thing_to_plot.msg_arrivals);
        title(sprintf([string_to_plot ': Message %d Path with Arrival Times'], msg_to_plot))
        if SAVE_BASE
            saveas(gcf, [test_str sprintf(['_msg_%d_arrivals_times_' string_to_plot '_' dataratesource '.png'], msg_to_plot)]); 
        end
    end
end


%% Functions
function dt = travel_time_table(j, k, i, A0, m_sizes)
    % returns travel time of message size m_i [m_sizes incorporated] along edge j->k
    % slant range - data rate table
    d = A0(j,k); % km
    if d == 0
        r = 0;
    elseif d <= 5
        r = 50; % kbps
    elseif d <= 300
        r = 10; % kbps
    elseif d <= 6000
        r = 1; % kbps
    else
        r = 0.125; % kbps
    end
    dt = m_sizes(i)/r; 
end
