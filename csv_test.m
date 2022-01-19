% This script is an example of how to use the function dataToGraph, which
% returns a container of adjacency matrices for the satellite constellation
% at each timestep. This script also plots each graph and converts the
% plots to an animated GIF.
% 
% Author: David Dezell Turner

%% Load variables

clear all
close all

load("ATrainCTrainGraphList");
load("GPSWalkerApproxGraphList");
load("MMSGraphList");
load("LunaNetGraphList");
load("LunaNetWithLLOGraphList");
load("DisconnectAttemptGraphList");

% graphList = LunaNetWithLLOGraphList;
% keyList = keys(graphList);
% filename = "LunaNetWithLLO_SecondAttempt_Animated.gif";

% accessFile = "DisconnectAttemptAccessData.xlsx";
% filename = strcat(erase(accessFile,"AccessData.xlsx"),'Animated.gif');
% graphList = dataToGraph(accessFile)
% keyList = keys(graphList);

%% Percent Time Connected Function test
at = percTimeConnected(ATrainCTrainGraphList)
gp = percTimeConnected(GPSWalkerApproxGraphList)
mm = percTimeConnected(MMSGraphList)
lnll = percTimeConnected(LunaNetWithLLOGraphList)
ln = percTimeConnected(LunaNetGraphList)
dis = percTimeConnected(DisconnectAttemptGraphList)

%% Shortest path tree test
% currentGraph = graph(graphList(keyList{9}));
% TR = shortestpathtree(currentGraph,1,'OutputForm','cell')
% tree = shortestpathtree(currentGraph,1);
% p = plot(currentGraph);
% highlight(p,tree,'EdgeColor','r');

%% Plot graphs and convert to GIF
% delay = 0.1;
% h = figure;
% axis tight manual
% 
% for i = 1:length(graphList)
%     ind = keyList{i};
%     currentGraph = graph(graphList(ind));
%     plot(currentGraph,"EdgeLabel",currentGraph.Edges.Weight)
% %     plot(currentGraph)
%     title(strcat("Time since epoch: ",num2str(ind)," s"))
%     drawnow
% 
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256);
% 
%     if i == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay); 
%     else 
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay); 
%     end 
% 
% end

%% Find min. distance (Optional, purely illustrative)
% minDist = Inf;
% minInd = 0;
% minMat = graphList(keyList{1});
% 
% for i = 1:length(graphList)
%     ind = keyList{i};
%     currentMat = graphList(ind);
%     if minDist > min(currentMat(currentMat>0))
%         minDist = min(currentMat(currentMat>0));
%         minMat = currentMat;
%     end
% end
% disp(minDist)
% disp(minInd)
% disp(minMat)