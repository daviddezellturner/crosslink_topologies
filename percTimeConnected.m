function percConnect = percTimeConnected(graphList)
    keyList = keys(graphList);
    percConnect = 0;
    for i = 1:length(graphList)
        ind = keyList{i};
        currentGraph = graph(graphList(ind));
        [bin,binsize] = conncomp(currentGraph);
        if binsize == numnodes(currentGraph)
            percConnect = percConnect + 1;
        end
    end
    percConnect = percConnect/length(graphList);
end