% Using STK, this code generates a constellation around a specified central
% body, then creates a CSV file of the distance between each pair of
% satellites at each timestep -- but only during time intervals where the
% pair of satellites have a direct line of sight between them.
%
% User defines constellation in "Define Constellation Parameters" section.
% 
% Author: David Dezell Turner

clear all
close all

%% Start STK
disp("Opening STK...")
app = actxserver('STK12.application');
app.UserControl = 1;

global root
root = app.Personality2;

global scenario
scenario = root.Children.New('eScenario','Lunar_Constellation_MATLAB_test');
scenario.SetTimePeriod('24 Feb 2012 12:00:00.000','25 Feb 2012 12:00:00.000');
scenario.StartTime = '24 Feb 2012 12:00:00.000';
scenario.StopTime = '25 Feb 2012 12:00:00.000';

root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec');

root.ExecuteCommand('Animate * Reset');
% How do I set the window to a Non-Earth-fixed view?

global satContainer
satContainer = containers.Map('KeyType','char','ValueType','any');

%% Define Constellation Parameters

% USER DEFINES CONSTELLATION HERE. Select a central body and call
% createConstellation() to define a constellation around it.

global satCentralBody;
satCentralBody = 'Moon'; % central body of satellite orbits

satsPerPlane = 3;
numPlanes = 1;
satName = "RingTest";
periAlt = 1000; % periapsis altitude [km]
apoAlt = 1000; % apoapsis altitude [km]
inc = 45; % [deg]
argPeri = 12; % argument of perigee [deg]
ascNode = 15; % RAAN
WalkerType = 'Delta';

createWalker(satName,numPlanes,satsPerPlane,periAlt,apoAlt,inc,argPeri,ascNode,WalkerType)
createWalker(satName,numPlanes,satsPerPlane,periAlt/2,apoAlt/2,inc,argPeri,ascNode+30,WalkerType)

%% Create Ground Stations
% One ground station for each DSN complex. Each is an approximation
% based on p. 15 of https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301K.pdf

% disp("Creating ground stations...")
% 
% dsn_gold = scenario.Children.NewOnCentralBody('eFacility','DSN_Goldstone','Earth');
% dsn_gold.Position.AssignGeodetic(35,243,1000); % lat, lon, alt [m]
% % dsn_gold.Graphics.Color = 3145645; % Yellow-Green
% 
% dsn_mad = scenario.Children.NewOnCentralBody('eFacility','DSN_Madrid','Earth');
% dsn_mad.Position.AssignGeodetic(40,355,800); % lat, lon, alt [m]
% % dsn_mad.Graphics.Color = 3145645; % Yellow-Green
% 
% dsn_can = scenario.Children.NewOnCentralBody('eFacility','DSN_Canberra','Earth');
% dsn_can.Position.AssignGeodetic(-35,148,700); % lat, lon, alt [m]
% % dsn_can.Graphics.Color = 3145645; % Yellow-Green
% 
% groundNameList = ["DSN_Goldstone","DSN_Madrid","DSN_Canberra"];
% groundStations = {dsn_gold,dsn_mad,dsn_can};
% groundContainer = containers.Map(groundNameList,groundStations);

%% Access Analysis

% accessContainer = containers.Map('KeyType','char','ValueType','any');
filename = strcat(satName,'AccessData.xlsx');
sheet = 1;

% Between satellites
allSats = keys(satContainer);
masterTimes = []; % Array of all timesteps for which there is access data
for a = 1:length(allSats)
    for b = a+1:length(allSats)
        try
            sat1 = satContainer(allSats{a});
            sat2 = satContainer(allSats{b});
    
            access = sat1.GetAccessToObject(sat2);
            access.ComputeAccess;
            
            accessName = strcat(allSats{a},"-to-",allSats{b});
    
            timeStep = 60;
            
            accessDP = access.DataProviders.Item('Access Data').Exec(scenario.StartTime,scenario.StopTime);
            accessStartTimes = cell2mat(accessDP.DataSets.GetDataSetByName('Start Time').GetValues);
            accessStopTimes = cell2mat(accessDP.DataSets.GetDataSetByName('Stop Time').GetValues);
            disp(strcat("Computing ",accessName,"..."))
    
            accessAER = access.DataProviders.Item('AER Data').Group.Item('BodyFixed').Exec(scenario.StartTime, scenario.StopTime, timeStep);
            AERTimes = cell2mat(accessAER.Interval.Item(cast(0, 'int32')).DataSets.GetDataSetByName('Time').GetValues);
            range = cell2mat(accessAER.Interval.Item(cast(0, 'int32')).DataSets.GetDataSetByName('Range').GetValues);
            for i = 1:1:accessAER.Interval.Count-1
                AERTimes = [AERTimes; cell2mat(accessAER.Interval.Item(cast(i, 'int32')).DataSets.GetDataSetByName('Time').GetValues)];
                range = [range; cell2mat(accessAER.Interval.Item(cast(i, 'int32')).DataSets.GetDataSetByName('Range').GetValues)];
            end
            AERTimes = timeStep*round(AERTimes./timeStep); % round time steps so they're all the same

            masterTimes = unique([masterTimes;AERTimes]);
    
    %             warning('off','MATLAB:xlswrite:AddSheet');
    %             xlswrite(filename,{access_name},sheet,'A1');
            headers = ["Time since start [s]","Distance between objects [km]"];
    %             xlswrite(filename,headers,sheet,'A2');
    %             xlswrite(filename,AERTimes,sheet,'A3');
    %             xlswrite(filename,range,sheet,'B3');
            writematrix(headers,filename,'Sheet',accessName,'Range','A1:B1');
            writematrix(AERTimes,filename,'Sheet',accessName,'Range','A2');
            writematrix(range,filename,'Sheet',accessName,'Range','B2');
%             sheet = sheet + 1;
        catch
            warning(strcat("No line of sight between ",allSats{a}," & ",allSats{b}));
        end
    end
end
writematrix(["Master Timestep List [s]"],filename,'Sheet','Master')
writematrix(masterTimes,filename,'Sheet','Master','Range','A2');
writematrix(["Master Satellite List"],filename,'Sheet','Master','Range','B1')
masterSatList = reshape(allSats,[length(allSats),1]);
writecell(masterSatList,filename,'Sheet','Master','Range','B2');

% % Between sats and ground stations
% for a = 1:length(allSats)
%     for b = 1:length(groundNameList)
%         try
%             sat1 = satContainer(allSats{a});
%             gs = groundContainer(groundNameList(b));
% 
%             access = sat1.GetAccessToObject(gs);
%             access.ComputeAccess;
%             accessDP = access.DataProviders.Item('Access Data').Exec(scenario.StartTime,scenario.StopTime);
%             accessStartTimes = accessDP.DataSets.GetDataSetByName('Start Time').GetValues;
%             accessStopTimes = accessDP.DataSets.GetDataSetByName('Stop Time').GetValues;
% 
%             access_name = strcat(allSats{a},"to",groundNameList(b));
%             disp(strcat("Computing ",access_name,"..."))
% 
%             xlswrite(filename,{access_name},sheet,'A1');
%             headers = {"Start Time","Stop Time"};
%             xlswrite(filename,headers,sheet,'A2');
%             xlswrite(filename,accperi_altessStartTimes,sheet,'A3');
%             xlswrite(filename,accessStartTimes,sheet,'B3');
%             sheet = sheet + 1;
%         catch
%             warning(strcat("No line of sight between ",allSats{a}," & ",groundNameList(b)));
%         end
%     end
% end

%% Functions
function createWalker(satName,numPlanes,satsPerPlane,periAlt,apoAlt,inc,argPeri,ascNode,WalkerType)
% CREATEWALKER Define a Walker constellation in STK and add the satellites
% to satContainer.
% 
% CREATEWALKER(satName,numPlanes,satsPerPlane,periAlt,apoAlt,inc,argPeri,ascNode,type):
%   satName = name of constellation [str]
%   numPlanes = number of planes [int]
%   satsPerPlane = number of satellites in each plane [int]
%   periAlt = altitude at periapsis, km [float]
%   apoAlt = altitude at apoapsis, km [float]
%   inc = inclination, deg [float]
%   argPeri = argument of periapsis, deg [float]
%   ascNode = right ascension of the ascending node, deg [float]
%   WalkerType = constellation type, "Delta" or "Star" [str]

    global satCentralBody;
    global satContainer;
    global root;
    global scenario;

    for a = 1:numPlanes
        for b = 1:satsPerPlane            
            satCurrentName = num2str(length(satContainer) + 1);
            disp(strcat("Creating Satellite ",satCurrentName,"..."))
    
            satellite = scenario.Children.NewOnCentralBody('eSatellite',satCurrentName,satCentralBody);
    
            keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
            keplerian.SizeShapeType = 'eSizeShapeAltitude';  % Changes from Ecc/Inc to Perigee/Apogee Altitude
            keplerian.LocationType = 'eLocationTrueAnomaly'; % Makes sure True Anomaly is being used
            keplerian.Orientation.AscNodeType = 'eAscNodeRAAN'; % RAAN (instead of LAN)
    
            keplerian.SizeShape.PerigeeAltitude = periAlt;
            keplerian.SizeShape.ApogeeAltitude = apoAlt;
            
            keplerian.Orientation.Inclination = inc;
            keplerian.Orientation.ArgOfPerigee = argPeri;
    
            if numPlanes == 1
                keplerian.Orientation.AscNode.Value = ascNode;
            else
                if WalkerType == "Delta"
                    keplerian.Orientation.AscNode.Value = mod(ascNode + 360*a/numPlanes, 360);
                elseif WalkerType == "Star"
                    keplerian.Orientation.AscNode.Value = mod(ascNode + 180*a/numPlanes, 360);
                else
                    error('WalkerType should be string with value "Delta" or "Star"')
                end
            end
            keplerian.Location.Value = 360*b/satsPerPlane;
    
            satellite.Propagator.InitialState.Representation.Assign(keplerian);
            satellite.Propagator.Propagate;
    
    %         satellite.Graphics.Attributes.Color = 16711680; % Blue
    
            satContainer(satCurrentName) = satellite;
    
%             % Generate ephemeris files for each satellite
%             saveFolder = pwd;
%             ephemFile = strcat(pwd,"\",satCurrentName,".e");
%             ephemText = 'ExportDataFile */Satellite/%s Ephemeris "%s" Type STK CoordSys J2000 CentralBody %s InterpBoundaries Include';
%             root.ExecuteCommand(sprintf(ephemText,satCurrentName,ephemFile,satCentralBody));
        end
    end
end