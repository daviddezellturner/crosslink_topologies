% This script is an example of how to use the function dataToGraph, which
% returns a container of adjacency matrices for the satellite constellation
% at each timestep. This script also plots each graph and converts the
% plots to an animated GIF.
% 
% Author: David Dezell Turner

clear all
close all

load("ATrainCTrainGraphList");
graphList = ATrainCTrainGraphList;
keyList = keys(graphList);

% graphList = dataToGraph("GPSWalkerApproxAccessData.xlsx")
% keyList = keys(graphList);

% Plots graphs and converts to GIF
delay = 0.1
h = figure;
axis tight manual
filename = 'ATrainCTrainAnimated.gif';

for i = 1:length(graphList)
    ind = keyList{i};
    currentGraph = graph(graphList(ind));
    plot(currentGraph,"EdgeLabel",currentGraph.Edges.Weight)
%     plot(currentGraph)
    title(strcat("Time since epoch: ",num2str(ind)," s"))
    drawnow

    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);

    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay); 
    end 

end