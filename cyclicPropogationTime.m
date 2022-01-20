function time = cyclicPropogationTime(graph)

    % (i don't think this is actually quite right - in reality the minimum
    % datrate will introduce a bottleneck. will still test with this in the
    % meantime..)

    datarates = distsToRates(graph)
    [greedy_cycle, raw_cost] = cyclefinder_backtrack_v2(graph,1);
    time = 0;
    % scale by node # because we are assuming aggregation of messages
    for node=1:length(greedy_cycle)-1
        time = time + node/datarates(greedy_cycle(node),greedy_cycle(node+1));
        time
    end
end