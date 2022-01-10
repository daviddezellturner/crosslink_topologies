function graphList = dataToGraph(filename, varargin)
% DATATOGRAPH Converts the constellation spreadsheet created by
% lunar_constellation_test.m into a container of adjacency matrices
% representing the constellation at each timestep.
% 
% DATATOGRAPH(filename):
%   filename = string, name of .xlsx file created by
%       lunar_constellation_test.m
%
% DATATOGRAPH(filename, numSteps):
%   filename = string, name of .xlsx file created by
%       lunar_constellation_test.m
%   numSteps = number of timesteps to use
%
% Author: David Dezell Turner

    sheets = sheetnames(filename);
    sheets = sheets(sheets ~= "Master");
    data = readmatrix(filename,'Sheet','Master','NumHeaderLines',1);
    timesteps = data(:,1);
    satNames = rmmissing(data(:,2));

    graphList = containers.Map('KeyType','double','ValueType','any');
    
    if nargin == 1
        final = length(timesteps);
    else
        final = varargin{1};
    end

    for i = 1:final
        disp(strcat("Extracting Timestep ",num2str(i)," out of ",num2str(final)))
        A = zeros(length(satNames)); % adjacency matrix

        for sheetnum = 1:length(sheets)
            satNums = split(sheets(sheetnum),"-to-");
            a = str2num(satNums(1));
            b = str2num(satNums(2));
            
            accessData = readmatrix(filename,'Sheet',sheets(sheetnum),'NumHeaderLines',1);
            accessTimes = accessData(:,1);
            accessDist = accessData(:,2);
            if any(accessTimes ==  timesteps(i))
                ind = find(accessTimes ==  timesteps(i));
                if length(ind) > 1
                    ind = ind(end);
                end
                A(a,b) = accessDist(ind);
                A(b,a) = accessDist(ind);
            end
        end
        graphList(timesteps(i)) = A;
    end
    disp("Adjacency matrices complete.")
end