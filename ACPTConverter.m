%%% ACPT Converter converts pre-sorted profile data into straight lines for
%%% the Alaska Coastal Profile Tool

%%% INPUTS
% csv with columns organized into 
% 01-05 Location|Year|Month|Day|Northing|
% 06-10 Easting|Datum|Elevation_Plot|Vertical_Datum_Plot|Line|
% 11-15 Point|Collector|Method|Grainsize|Grainsize_Method|
% 16-20 Vertical_Uncertainty|Horizontal_Uncertainty|Elevation_Raw|Vertical_Datum_Raw|Notes
% nodata values:
%   str = none
%   dbl = -999
% no empty cells

%%% OUTPUT
% One csv for each Line value, 
% sorted by detected seaward-to-landward point number,
% straightened
% Uncertainty csv with columns for
% Profile number
% max offset (meters)
% mean offset (meters)
% standard deviation (meters)

%%% Written by Richard Buzard, 9/3/2018

%% Setup
clear all; close all; clc
addpath(cd)
strings = [1,7,9,12,13,15,19,20];   % variable used much later
[CSVfile, CSVpath] = uigetfile({'*.csv'}, 'Select CSV', 'Multiselect', 'On');
cd(CSVpath)
%newPath = uigetdir(CSVpath,'Select Output Location');
CSV_separate = [];

% Delete files if they exist
try cd(strcat([fileparts(pwd), '\SEPARATE']))
    oldFiles = dir(cd);
    for ff = 1:length(oldFiles)
        if strcmp(oldFiles(ff).name(5:10), CSVfile(1:6))
            delete(oldFiles(ff).name)
        end
    end
catch ME   
end
try cd(strcat([fileparts(pwd), '\STRAIGHT']))
    oldFiles = dir(cd);
    for ff = 1:length(oldFiles)
        if strcmp(oldFiles(ff).name(5:10), CSVfile(1:6))
            delete(oldFiles(ff).name)
        end
    end
catch ME   
end
disp('Working...')
%% Import
CSVfile = cellstr(CSVfile);
for kk = 1: numel(CSVfile)          % For each file
    cd(CSVpath)    
    try fid = fopen(CSVfile{kk});       % Make file id
    catch ME
        fid = fopen(CSVfile);
    end
% Import the datasheet as a cell array
c = textscan(fid, '%s %f %f %f %f %f %s %f %s %f %f %s %s %f %s %f %f %f %s %s','Delimiter',',');

% Detect and remove headers
if isempty(c{2})
    c = textscan(fid, '%s %f %f %f %f %f %s %f %s %f %f %s %s %f %s %f %f %f %s %s','Delimiter',',','HeaderLines',1);
end

% Detect incomplete datasheets
for hh = 1:20
    if isempty(c{hh})
        error(strcat(['Error: Column ',num2str(hh), ' is wrong format or empty.']))
    end
end

fclose(fid);                    % Close file

% work in the placement folder
try cd(strcat([fileparts(pwd), '\SEPARATE']))
catch ME
    cd(fileparts(pwd))
    mkdir('SEPARATE')
    cd(strcat([pwd, '\SEPARATE']))
    
end
pathSeparate = cd;

m = [];                         % Reset m matrix

% Change Grain Size nodata values to -999
for gg = 1:length(c{14})
    if c{14}(gg) == 0
       c{14}(gg) = -999;
    end
    if isempty(char(c{20}(gg)))
       c{20}(gg) = {'none'};
    end
end

% Slice Those Lines

for ii = 0:max(c{10})           % for each line value        
    for jj = 1:length(c{10})    % for each row of the datasheet
        if ii == c{10}(jj)      % if the value of "line" is this line
                                % put the data into a matrix
           p = cellfun(@(v) v(jj), c(1,:),'UniformOutput',0);
           for ss = strings     % remove embedded cell arrays
                p{1,ss} = char(p{1,ss});
           end
           m = [m;p];           % tack on the line
        end
    end
    
    if isempty(m) ~= 1
        
        % Check to make sure transect goes SL to land
        fitvars = polyfit([m{:,11}], [m{:,8}], 1);
        if  fitvars(1) < -0  
            disp(strcat(num2str(c{2}(1)*1000+c{3}(1)*100+c{4}(1)),...
                {' Check that Profile '}, num2str(ii),{' direction is sea to land'}))
            % this used to automatically flip the Z values assuming the
            % data were input wrong. However it is better that the user
            % verify the data input is correct.
%             n = length([m{:,11}]);
%             for aa = 1:n
%                 m{aa,11} = n + 1 - m{aa,11};
%             end
        end

        m = sortrows(m,11);             % Sort the rows by point
        
        % Check to make sure transect goes in order
        x = [m{:,6}]; xOffset = x(1);
        y = [m{:,5}]; yOffset = y(1);
        x1 = x - xOffset;
        y1 = y - yOffset;
        xypos = sqrt(x1.^2 + y1.^2); 
        xydif = xypos(2:end)-xypos(1:end-1);
        if min(xydif)<-2        % if there are out-of-order points further than 2 meters
            [~,I] = sort(xypos);            % sort by along-profile distance
            m = m(I,:);
            for aa = 1:size(m,1)
                m{aa,11} = aa;
            end
        end
    

        % Put the matrix into a new csv
        year = num2str(m{1,2}, '%04i');
        month = num2str(m{1,3}, '%02i');
        day = num2str(m{1,4}, '%02i');
        if or(strcmp(m{1,13},'DGPS'),strcmp(m{1,13},'GPS'))
             method = 'GPS';
        else method = 'RMS';
        end
        lineNum = num2str(ii, '%03i');
        filename = strcat([lineNum, '_', year, month, day, '000000', method, '.csv']);
        CSV_separate = [CSV_separate; filename]; 
        cell2csv(filename, m, ',', '2013', '.')
    end
    % close it down and do it again
    m = [];
    close all
end
end
%% Profile Solver Alpha

%addpath(cd)
%[CSVfile, CSVpath] = uigetfile({'*.csv'}, 'Select CSV', 'Multiselect', 'On');
CSVfile = cellstr(CSV_separate); % make cell array of CSVfile(s) in case there is just one
numProf = numel(CSVfile);   % count number of files
[data, coeff, D, unc, xc, yc, z, pos] = deal(1:numProf) ;   % set up variables for each file
strings = [1,7,9,12,13,15,19,20];   % variable used much later

%% run solver
for kk = 1: numProf        % For each file
    try fid = fopen(strcat(pathSeparate,'\',CSVfile{kk}));       % Make file id
    catch ME
        fid = fopen(CSVfile);
    end
    % Import the datasheet as a cell array
    c = textscan(fid, '%s %f %f %f %f %f %s %f %s %f %f %s %s %f %s %f %f %f %s %s','Delimiter',',');
    fclose(fid);                    % Close file
    
    % Work in placement folder
    try cd(strcat([fileparts(pwd), '\STRAIGHT']))
        catch ME
    cd(fileparts(pwd))
    mkdir('STRAIGHT')
    cd(strcat([pwd, '\STRAIGHT']))
    end

    

% Find or make directory to store corrected profiles
% % % try cd(strcat([fileparts(pwd), '\ProfilesCorrected']))
% % % catch ME
% % %     cd(fileparts(pwd))
% % %     mkdir('ProfileUncertainty')
% % %     cd(strcat([pwd, '\ProfileUncertainty']))
% % % end
    

% Solve the Profile for linear fit
format long

    % Make numbers smaller so Matlab can calculate stuff
    xOffset = round(min(c{1,6}),-4);
    yOffset = round(min(c{1,5}),-4);

    x = c{1,6} - xOffset;
    y = c{1,5} - yOffset;
    z = c{1,8};
    [coeff,S] = polyfit(x,y,1);
    
    % Find distance perpendicular to linear fit
    A = coeff(1);
    B = -1;
    C = coeff(2);
    D = (A*x + B*y + C)/sqrt(A^2 + B^2);
    
    % Solve for corrected position
    theta = atan(A/B);
    xd = D*sin(theta);
    yd = D*cos(theta);
    
    xc = (x + xd) + xOffset;
    yc = (y + yd) + yOffset; 
    
    % Setup matrix with corrected data
    c{1,6} = xc;
    c{1,5} = yc;
    
    % Clean matrix for saving
    m =[];
    for jj = 1:length(c{10})    % for each row of the datasheet
                                % put the data into a matrix
           p = cellfun(@(v) v(jj), c(1,:),'UniformOutput',0);
           for ss = strings     % remove embedded cell arrays
                p{1,ss} = char(p{1,ss});
           end
           m = [m;p];           % tack on the line
    end

    
    
    % Save results as new .csv
    filename = char(CSVfile{kk});
    cell2csv(filename, m, ',', '2013', '.')
    c = {};
    close all
    % Make log file for uncertainty values
     
    Uncert(kk,1) = str2double(CSVfile{kk}(1:3));
    Uncert(kk,2) = max(abs(D));
    Uncert(kk,3) = mean(D);
    Uncert(kk,4) = std(D);
end
%% Store uncertainty
% Find or make directory to store uncertainties
try cd(strcat([fileparts(pwd), '\UNCERTAINTY']))
catch ME
    cd(fileparts(pwd))
    mkdir('UNCERTAINTY')
    cd(strcat([pwd, '\UNCERTAINTY']))
end

% Create uncertainty CSV
csvwrite(strcat(['Uncert',CSVfile{1,1}(4:end)]),Uncert)
fclose('all');
cd(CSVpath)
disp('Success!')


%% Internal Functions
function cell2csv(fileName, cellArray, separator, excelYear, decimal)
% Writes cell array content into a *.csv file.
% 
% CELL2CSV(fileName, cellArray, separator, excelYear, decimal)
%
% fileName     = Name of the file to save. [ i.e. 'text.csv' ]
% cellArray    = Name of the Cell Array where the data is in
% separator    = sign separating the values (default = ';')
% excelYear    = depending on the Excel version, the cells are put into
%                quotes before they are written to the file. The separator
%                is set to semicolon (;)
% decimal      = defines the decimal separator (default = '.')
%
%         by Sylvain Fiedler, KA, 2004
% updated by Sylvain Fiedler, Metz, 06
% fixed the logical-bug, Kaiserslautern, 06/2008, S.Fiedler
% added the choice of decimal separator, 11/2010, S.Fiedler

%% Checking für optional Variables
if ~exist('separator', 'var')
    separator = ',';
end

if ~exist('excelYear', 'var')
    excelYear = 1997;
end

if ~exist('decimal', 'var')
    decimal = '.';
end

%% Setting separator for newer excelYears
if excelYear > 2000
    separator = ';';
end

%% Write file

datei = fopen(fileName, 'w');

for z=1:size(cellArray, 1)
    for s=1:size(cellArray, 2)
        
        var = eval(['cellArray{z,s}']);
        % If zero, then empty cell
        if size(var, 1) == 0
            var = '';
        end
        % If numeric -> String
        if isnumeric(var)
            var = num2str(var);
            % Conversion of decimal separator (4 Europe & South America)
            % http://commons.wikimedia.org/wiki/File:DecimalSeparator.svg
            if decimal ~= '.'
                var = strrep(var, '.', decimal);
            end
        end
        % If logical -> 'true' or 'false'
        if islogical(var)
            if var == 1
                var = 'TRUE';
            else
                var = 'FALSE';
            end
        end
        % If newer version of Excel -> Quotes 4 Strings
        if excelYear > 2000
            var = ['"' var '"'];
        end
        
        % OUTPUT value
        fprintf(datei, '%s', var);
        
        % OUTPUT separator
        if s ~= size(cellArray, 2)
            fprintf(datei, separator);
        end
    end
    if z ~= size(cellArray, 1) % prevent a empty line at EOF
        % OUTPUT newline
        fprintf(datei, '\n');
    end
end
% Closing file
fclose(datei);
end