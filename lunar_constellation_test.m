% Using STK, this code generates a constellation around a specified central
% body, then creates a CSV file of the distance between each pair of
% satellites at each timestep -- but only during time intervals where the
% pair of satellites have a direct line of sight between them.
%
% User defines constellation in "Define Constellation Parameters" section.
% 
% Author: David Dezell Turner


% For setting colors: https://docs.microsoft.com/en-us/office/vba/api/excel.xlrgbcolor

clear all
close all

%% Start STK
disp("Opening STK...")
app = actxserver('STK12.application');
app.UserControl = 1;

global root
root = app.Personality2;

global scenario

start = '24 Feb 2022 12:00:00.000';
stop = '25 Feb 2022 12:00:00.000';
satName = "LunaNet";
timeStep = 60; %43200; % [s]

scenario = root.Children.New('eScenario','Lunar_Constellation_MATLAB_test');
scenario.SetTimePeriod(start,stop);
scenario.StartTime = start;
scenario.StopTime = stop;

root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec');

root.ExecuteCommand('Animate * Reset');
% How do I set the window to a Non-Earth-fixed view?

global satContainer
satContainer = containers.Map('KeyType','char','ValueType','any');

%% Define Constellation Parameters

% USER DEFINES CONSTELLATION HERE. Select a central body and call
% createConstellation() to define a constellation around it.

global satCentralBody;
satCentralBody = 'Earth'; % central body of satellite orbits

satsPerPlane = 4;
numPlanes = 6;
periAlt = 20200; % periapsis altitude [km]
apoAlt = 20200; % apoapsis altitude [km]
inc = 55; % [deg]
argPeri = 30.442; % argument of perigee [deg]
ascNode = 14.28; % RAAN
WalkerType = 'Delta';

MMS = ["40482","40483","40484","40485"]; % Magnetospheric Multiscale
ATrain = ["27424","28376","38337","40059"]; % Earth Observing System (EOS) A-Train
CTrain = ["29108","29107"];
GPS = ["39533","39741","40105","40294","40534","40730","41019","41328","43873","44506","45854","46826","48859","24876","27704","28129","32711","35752","26360","29486","28874","32260","27663","32384","29601","28190","28474","36585","37753","38833","39166"];

createLunaNet()
% satFromDB(MMS)
% createWalker(numPlanes,satsPerPlane,periAlt,apoAlt,inc,argPeri,ascNode,WalkerType)
% createWalker(numPlanes,satsPerPlane,periAlt/2,apoAlt/2,inc,argPeri,ascNode+30,WalkerType)

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

function satFromDB(SSCNumList)
    global satCentralBody;
    global satContainer;
    global root;
    global scenario;

    for i = 1:length(SSCNumList)
        satCurrentName = num2str(length(satContainer) + 1);
        disp(satCurrentName)
        try
            root.ExecuteCommand(strcat("ImportFromDB * Satellite AGIServer Propagate On TimePeriod UseScenarioInterval SSCNumber ",SSCNumList(i), " Rename ", satCurrentName));
            satContainer(satCurrentName) = scenario.Children.Item(int32(i-1));
        catch
            warning(strcat(SSCNumList(i)," not imported."))
        end
    end
end

function createWalker(numPlanes,satsPerPlane,periAlt,apoAlt,inc,argPeri,ascNode,WalkerType)
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

function createLunaNet()
    global satCentralBody;
    global satContainer;
    global root;
    global scenario;

    satCentralBody = 'Moon';
    semimajor = 6142.4; % [km]

    % 12-Hr Circular Equatorial Sat
    disp("Creating 12-Hr Circular Equatorial Sat...")
    satellite = scenario.Children.NewOnCentralBody('eSatellite','1',satCentralBody);
    keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
    keplerian.SizeShapeType = 'eSizeShapeSemimajorAxis';
    keplerian.LocationType = 'eLocationMeanAnomaly';
    keplerian.Orientation.AscNodeType = 'eAscNodeRAAN';
    keplerian.SizeShape.SemimajorAxis = semimajor;
    keplerian.SizeShape.Eccentricity = 0;
    keplerian.Orientation.Inclination = 0;
    keplerian.Orientation.AscNode.Value = 0;
    keplerian.Orientation.ArgOfPerigee = 315;
    keplerian.Location.Value = 0;
    graphics = satellite.Graphics;
    graphics.SetAttributesType('eAttributesBasic');
    attributes = graphics.Attributes;
    attributes.Color = 65535; % Yellow
    satellite.Propagator.InitialState.Representation.Assign(keplerian);
    satellite.Propagator.Propagate;
    satContainer("1") = satellite;
    
    ecc = 0.59999;
    inc = 57.7; % deg

    % 12-Hr Elliptical North Sat 1
    disp("Creating 12-Hr Elliptical North Sat 1...")
    satellite = scenario.Children.NewOnCentralBody('eSatellite',"2",satCentralBody);
    keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
    keplerian.SizeShapeType = 'eSizeShapeSemimajorAxis';
    keplerian.LocationType = 'eLocationTrueAnomaly';
    keplerian.Orientation.AscNodeType = 'eAscNodeRAAN';
    keplerian.SizeShape.SemimajorAxis = semimajor;
    keplerian.SizeShape.Eccentricity = ecc;
    keplerian.Orientation.Inclination = inc;
    keplerian.Orientation.AscNode.Value = 270;
    keplerian.Orientation.ArgOfPerigee = 270;
    keplerian.Location.Value = 0;
    graphics = satellite.Graphics;
    graphics.SetAttributesType('eAttributesBasic');
    attributes = graphics.Attributes;
    attributes.Color = 255; % Red
    satellite.Propagator.InitialState.Representation.Assign(keplerian);
    satellite.Propagator.Propagate;
    satContainer("2") = satellite;

    % 12-Hr Elliptical North Sat 2
    disp("Creating 12-Hr Elliptical North Sat 2...")
    satellite = scenario.Children.NewOnCentralBody('eSatellite',"3",satCentralBody);
    keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
    keplerian.SizeShapeType = 'eSizeShapeSemimajorAxis';
    keplerian.LocationType = 'eLocationTrueAnomaly';
    keplerian.Orientation.AscNodeType = 'eAscNodeRAAN';
    keplerian.SizeShape.SemimajorAxis = semimajor;
    keplerian.SizeShape.Eccentricity = ecc;
    keplerian.Orientation.Inclination = inc;
    keplerian.Orientation.AscNode.Value = 270;
    keplerian.Orientation.ArgOfPerigee = 270;
    keplerian.Location.Value = 180;
    graphics = satellite.Graphics;
    graphics.SetAttributesType('eAttributesBasic');
    attributes = graphics.Attributes;
    attributes.Color = 255; % Red
    satellite.Propagator.InitialState.Representation.Assign(keplerian);
    satellite.Propagator.Propagate;
    satContainer("3") = satellite;

    % 12-Hr Elliptical South Sat 1
    disp("Creating 12-Hr Elliptical South Sat 1...")
    satellite = scenario.Children.NewOnCentralBody('eSatellite',"4",satCentralBody);
    keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
    keplerian.SizeShapeType = 'eSizeShapeSemimajorAxis';
    keplerian.LocationType = 'eLocationTrueAnomaly';
    keplerian.Orientation.AscNodeType = 'eAscNodeRAAN';
    keplerian.SizeShape.SemimajorAxis = semimajor;
    keplerian.SizeShape.Eccentricity = ecc;
    keplerian.Orientation.Inclination = inc;
    keplerian.Orientation.AscNode.Value = 0;
    keplerian.Orientation.ArgOfPerigee = 90;
    keplerian.Location.Value = 0;
    graphics = satellite.Graphics;
    graphics.SetAttributesType('eAttributesBasic');
    attributes = graphics.Attributes;
    attributes.Color = 16711680; % Blue
    satellite.Propagator.InitialState.Representation.Assign(keplerian);
    satellite.Propagator.Propagate;
    satContainer("4") = satellite;

    % 12-Hr Elliptical South Sat 2
    disp("Creating 12-Hr Elliptical South Sat 2...")
    satellite = scenario.Children.NewOnCentralBody('eSatellite',"5",satCentralBody);
    keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
    keplerian.SizeShapeType = 'eSizeShapeSemimajorAxis';
    keplerian.LocationType = 'eLocationTrueAnomaly';
    keplerian.Orientation.AscNodeType = 'eAscNodeRAAN';
    keplerian.SizeShape.SemimajorAxis = semimajor;
    keplerian.SizeShape.Eccentricity = ecc;
    keplerian.Orientation.Inclination = inc;
    keplerian.Orientation.AscNode.Value = 0;
    keplerian.Orientation.ArgOfPerigee = 90;
    keplerian.Location.Value = 180;
    graphics = satellite.Graphics;
    graphics.SetAttributesType('eAttributesBasic');
    attributes = graphics.Attributes;
    attributes.Color = 16711680; % Blue
    satellite.Propagator.InitialState.Representation.Assign(keplerian);
    satellite.Propagator.Propagate;
    satContainer("5") = satellite;

end