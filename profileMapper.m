%%% Profile Mapper
% 
% This script creates a CSV of profile locations to use in Arc
% 
% INPUTS
% Straightened profiles produced by ACPTConverter
% 
% OUTPUTS
% CSV of profile locations to be put converted to lines in ArcGIS using XY
% to Line (geoprocessing tool in Arc)
% CSV also shows profile count and total years between collection
% 
% LIMITS
% The code is designed to take in unlimited repeat profiles, but it always
% compares them to the first profile collected
% 
% Written by Richard Buzard, February 1, 2022
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Setup
clear all; close all;

% select files
[CSVfile, CSVpath] = uigetfile({'*.csv'}, 'Select CSV files of straightened profiles', 'Multiselect', 'On');
CSVfile = cellstr(CSVfile);
numProf = numel(CSVfile);
strings = [1,7,9,12,13,15,19,20];   % variable used much later

%% extract xyz data
for kk = 1: numProf        % For each file
    cd(CSVpath)
    try fid = fopen(CSVfile{kk});       % Make file id
    catch ME
        fid = fopen(CSVfile);
    end
    % Import the datasheet as a cell array
    c{kk} = textscan(fid, '%s %f %f %f %f %f %s %f %s %f %f %s %s %f %s %f %f %f %s %s','Delimiter',',');
    fclose(fid);       
    pNums(kk) = c{kk}{1,10}(1);
    x{kk} = c{kk}{6};
    y{kk} = c{kk}{5};
    z{kk} = c{kk}{8};
    d(kk) = str2double(CSVfile{kk}(5:12));
end
[profNums, ia, ic] = unique(pNums);  % index unique profile numbers

%% Summarize profiles

% collapse profiles onto one line and store in variable
m=[];
for pp = 1:length(profNums)       % for each profile number
    ppidx = find(pNums==profNums(pp));          % index where the profiles exist
    pCount = numel(ppidx);                      % find how many profiles exist for it
    m(pp,1) = profNums(pp);                     % record profile number
    m(pp,6) = pCount;
    m(pp,7) = round(max(d(ppidx))/10000) - round(min(d(ppidx))/10000);
    startx=[]; starty=[]; endx=[]; endy=[];
    for kk = 1:pCount           % for each profile of this number
        x2 = x{ppidx(kk)};      % get xy 
        y2 = y{ppidx(kk)};
             % for each profile
        startx(kk) = x2(1); % find start points of profiles
        starty(kk) = y2(1);
        endx(kk)   = x2(end);
        endy(kk)   = y2(end);
        dirx = sign(x2(end) - startx(kk));  % get line direction
        diry = sign(y2(end) - starty(kk));
    end
    % determine max profile length coordinates
    if dirx == 1                % if profile goes east
        m(pp,2) = min(startx);  % start west
        m(pp,4) = max(endx);  % end east
    else                        % if profile goes west
        m(pp,2) = max(startx);  % start east
        m(pp,4) = min(endx);  % end west
    end
    if diry == 1                % if profile goes north
        m(pp,3) = min(starty);    % start south
        m(pp,5) = max(endy);    % end north
    else                        % if profile goes south
        m(pp,3) = max(starty);    % start north
        m(pp,5) = min(endy);    % end south
    end
   
end

%% Export as csv
t = array2table(m);
t.Properties.VariableNames(1:7) = {'ProfNum','startX','startY','endX','endY','ProfCount','YearsApart'};

try outfile = strcat(CSVpath(1:strfind(CSVpath,'ACPT')+4),'profileLocations.csv');
    disp(strcat({'File output: '},outfile))
catch ME
    [outfile, outpath] = uiputfile('*.csv','Choose output location and filename','profileLocations');
    outfile = fullfile(outpath,outfile);
end
writetable(t,outfile)
disp('In ArcGIS, use tool XY to Line to display profile locations')