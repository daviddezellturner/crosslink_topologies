% Author: David Dezell Turner
% 
% Using STK, this code generates a constellation around a specified central
% body, then creates a CSV file listing the time intervals during which
% each satellite or ground station is visible from each satellite.

clear all
close all

% Helpful resource: https://help.agi.com/stkdevkit/11.4.0/Content/stkObjects/ObjModMatlabCodeSamples.htm#111

%% Define Constellation Parameters
satCentralBody = 'Moon'; % central body of satellite orbits
satsPerPlane = 2;
numPlanes = 6;
satName = "LunarSat";
peri_alt = 1000; % periapsis altitude [km]
apo_alt = 1000; % apoapsis altitude [km]
inclination = 45; % [deg]
arg_perigee = 12; % argument of perigee [deg]

%% Start STK
disp("Opening STK...")
app = actxserver('STK12.application');
app.UserControl = 1;
root = app.Personality2;

scenario = root.Children.New('eScenario','Lunar_Constellation_MATLAB_test');
scenario.SetTimePeriod('24 Feb 2012 12:00:00.000','25 Feb 2012 12:00:00.000');
scenario.StartTime = '24 Feb 2012 12:00:00.000';
scenario.StopTime = '25 Feb 2012 12:00:00.000';

root.ExecuteCommand('Animate * Reset');
% How do I set the window to a Non-Earth-fixed view?

%% Create constellation
satContainer = containers.Map('KeyType','char','ValueType','any');

for a = 1:numPlanes
    for b = 1:satsPerPlane
        satCurrentName = strcat(satName,num2str(a),num2str(b));
        disp(strcat("Creating ",satCurrentName,"..."))
        satNameList(a,b) = satCurrentName;

        satellite = scenario.Children.NewOnCentralBody('eSatellite',satCurrentName,satCentralBody);

        keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
        keplerian.SizeShapeType = 'eSizeShapeAltitude';  % Changes from Ecc/Inc to Perigee/Apogee Altitude
        keplerian.LocationType = 'eLocationTrueAnomaly'; % Makes sure True Anomaly is being used
        keplerian.Orientation.AscNodeType = 'eAscNodeLAN'; % Use LAN instead of RAAN for data entry

        keplerian.SizeShape.PerigeeAltitude = peri_alt;      % km
        keplerian.SizeShape.ApogeeAltitude = apo_alt;       % km
        
        keplerian.Orientation.Inclination = inclination;         % deg
        keplerian.Orientation.ArgOfPerigee = arg_perigee;        % deg
        keplerian.Orientation.AscNode.Value = 360*a/numPlanes;       % deg
        keplerian.Location.Value = 360*b/satsPerPlane;                 % deg
        
        satellite.Propagator.InitialState.Representation.Assign(keplerian);
        satellite.Propagator.Propagate;

        satContainer(satCurrentName) = satellite;
    end
end

% root.ExecuteCommand('Walker */Satellite/LunarSat Type Delta NumPlanes 7 NumSatsPerPlane 2 InterPlanePhaseIncrement 1 ColorByPlane Yes')

%% Create Ground Stations
% One ground station for each DSN complex. Each is an approximation
% based on p. 15 of https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301K.pdf
disp("Creating ground stations...")

dsn_gold = scenario.Children.NewOnCentralBody('eFacility','DSN_Goldstone','Earth');
dsn_gold.Position.AssignGeodetic(35,243,1000); % lat, lon, alt [m]

dsn_mad = scenario.Children.NewOnCentralBody('eFacility','DSN_Madrid','Earth');
dsn_mad.Position.AssignGeodetic(40,355,800); % lat, lon, alt [m]

dsn_can = scenario.Children.NewOnCentralBody('eFacility','DSN_Canberra','Earth');
dsn_can.Position.AssignGeodetic(-35,148,700); % lat, lon, alt [m]

groundNameList = ["DSN_Goldstone","DSN_Madrid","DSN_Canberra"];
groundStations = {dsn_gold,dsn_mad,dsn_can};
groundContainer = containers.Map(groundNameList,groundStations);

%% Access Analysis

% accessContainer = containers.Map('KeyType','char','ValueType','any');
filename = strcat(satName,'AccessData.xlsx');
sheet = 1;

% Between satellites
allSats = keys(satContainer);
for a = 1:length(allSats)
    for b = a+1:length(allSats)
        try
            sat1 = satContainer(allSats{a});
            sat2 = satContainer(allSats{b});

            access = sat1.GetAccessToObject(sat2);
            access.ComputeAccess;
            accessDP = access.DataProviders.Item('Access Data').Exec(scenario.StartTime,scenario.StopTime);
            accessStartTimes = accessDP.DataSets.GetDataSetByName('Start Time').GetValues;
            accessStopTimes = accessDP.DataSets.GetDataSetByName('Stop Time').GetValues;

            access_name = strcat(allSats{a},"to",allSats{b});
            disp(strcat("Computing ",access_name,"..."))

            xlswrite(filename,{access_name},sheet,'A1');
            headers = {"Start Time","Stop Time"};
            xlswrite(filename,headers,sheet,'A2');
            xlswrite(filename,accessStartTimes,sheet,'A3');
            xlswrite(filename,accessStartTimes,sheet,'B3');
            sheet = sheet + 1;
        catch
            warning(strcat("Error: ",allSats{a}," & ",allSats{b}));
        end
    end
end

% Between sats and ground stations
for a = 1:length(allSats)
    for b = 1:length(groundNameList)
        try
            sat1 = satContainer(allSats{a});
            gs = groundContainer(groundNameList(b));

            access = sat1.GetAccessToObject(gs);
            access.ComputeAccess;
            accessDP = access.DataProviders.Item('Access Data').Exec(scenario.StartTime,scenario.StopTime);
            accessStartTimes = accessDP.DataSets.GetDataSetByName('Start Time').GetValues;
            accessStopTimes = accessDP.DataSets.GetDataSetByName('Stop Time').GetValues;

            access_name = strcat(allSats{a},"to",groundNameList(b));
            disp(strcat("Computing ",access_name,"..."))

            xlswrite(filename,{access_name},sheet,'A1');
            headers = {"Start Time","Stop Time"};
            xlswrite(filename,headers,sheet,'A2');
            xlswrite(filename,accessStartTimes,sheet,'A3');
            xlswrite(filename,accessStartTimes,sheet,'B3');
            sheet = sheet + 1;
        catch
            warning(strcat("Error: ",allSats{a}," & ",groundNameList(b)));
        end
    end
end