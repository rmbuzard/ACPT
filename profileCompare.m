%%% Profile Compare
% 
% This script plots coastal profile data
% 
% INPUTS
% Straightened profiles produced by ACPTConverter
% 
% OUTPUTS
% Interactive plot to examine profile data
% 
% LIMITS
% The code is designed to take in unlimited repeat profiles, but it always
% compares them to the first profile collected
% 
% Written by Richard Buzard, October 11, 2017
% Updated December 8, 2021
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Setup
clear all; close all;

% select files
[CSVfile, CSVpath] = uigetfile({'*.csv'}, 'Select CSV files of all profiles', 'Multiselect', 'On');
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
    %profNum(kk) = c{kk}{10}(1); % get profile number for this file
    x{kk} = c{kk}{6};
    y{kk} = c{kk}{5};
    z{kk} = c{kk}{8};
    d(kk) = str2double(CSVfile{kk}(5:12));
end

% Compare profiles
[profNums, ia, ic] = unique(pNums);  % index unique profile numbers

for pp = 1:length(profNums)       % for each profile number
    ppidx = find(pNums==profNums(pp));       % index where the profiles exist
    pCount = numel(ppidx);          % find how many profiles exist for it
    if pCount > 1                   % if there are >1 profiles to compare
        xOffset = round(min(x{ppidx(1)}),-4);   % get offset for polyfit calculation
        yOffset = round(min(y{ppidx(1)}),-4);
        x1 = x{ppidx(1)} - xOffset;     % get master profile xyz
        y1 = y{ppidx(1)} - yOffset;
        z1 = z{ppidx(1)};
        [coeff,S] = polyfit(x1,y1,1);   % find slope of first profile
        
        for jj = 2:pCount           % for each subsequent profile of this number
            x2 = x{ppidx(jj)} - xOffset;      % get xyz 
            y2 = y{ppidx(jj)} - yOffset;
            z2 = z{ppidx(jj)};
            A = coeff(1);           % find perpendicular distance D between profiles
            B = -1;
            C = coeff(2);
            D = (A*x2 + B*y2 + C)/sqrt(A^2 + B^2);
            theta = atan(A/B);      % Solve for corrected position
            xd = D*sin(theta);
            yd = D*cos(theta);
            xc = (x2 + xd) + xOffset; % put it back into the coordinate system
            yc = (y2 + yd) + yOffset; 
            profiles{pp}{jj}.date = d(ppidx(jj));   % store corrections
            profiles{pp}{jj}.x    = xc;
            profiles{pp}{jj}.y    = yc;
            profiles{pp}{jj}.z    = z2;
        end
    end
    profiles{pp}{1}.date = d(ppidx(1));
    profiles{pp}{1}.x    = x{ppidx(1)};
    profiles{pp}{1}.y    = y{ppidx(1)};
    profiles{pp}{1}.z    = z{ppidx(1)};
end

% calculate elevation changes

for pp = 1:length(profiles)             % for each profile
    for kk = 1:length(profiles{pp})         % for each profile dataset
        startx(kk) = profiles{pp}{kk}.x(1); % find start points of profiles
        starty(kk) = profiles{pp}{kk}.y(1);
        dirx = sign(profiles{pp}{kk}.x(end) - startx(kk));  % get line direction
        diry = sign(profiles{pp}{kk}.y(end) - starty(kk));
    end
    % determine which profile starts furthest seaward
    if length(startx)>1 % if there is more than one acquisition of this profile
        if sign(startx(2) - startx(1)) == dirx  % find profile start point furthest seaward
            xOffset = startx(1);
        else
            xOffset = startx(2);
        end
        if sign(starty(2) - starty(1)) == diry
            yOffset = starty(1);
        else
            yOffset = starty(2);
        end
    end
    % get elevations along profile from start point
    
    for kk = 1:length(profiles{pp})
        x1 = profiles{pp}{kk}.x - xOffset;
        y1 = profiles{pp}{kk}.y - yOffset;
        profiles{pp}{kk}.pos = sqrt(x1.^2 + y1.^2);
    end
   
    
%    for kk = 1:length(profiles{pp})     % check the plot to see if they line up
%         plot(profiles{pp}{kk}.x,profiles{pp}{kk}.y)
%         hold on
%    end
end
%% Plot all profiles in big subplot
clf
for pp = 1:length(profiles)
    subplot(round(length(ia)/2),2,pp)
    for kk = 1:length(profiles{pp})
        plot(profiles{pp}{kk}.pos,profiles{pp}{kk}.z)
        hold on
        ax = gca;
        if or(pp == length(ia), pp == length(ia)-1)
        else
            ax.XTickLabel = [];
        end
    end
    title(profNums(pp))
    ylabel(profNums(pp),'rotation',0)
    axis equal
end



%% Plot interactive graph
itemNames ={};
for pp = 1:length(profiles)
    n = strcat({'Profile '},num2str(profNums(pp)),...
        {' ('},num2str(length(profiles{pp})),{' dates)'});
    itemNames{pp,1} = n{1};
end
pp = 1;

plotOptions(profNums,profiles,pp,itemNames);

%% FUNCTIONS
function plotOptions(profNums,profiles,pp,itemNames)
fig = uifigure;
fig.Position(3:4) = [430 330];
ax = uiaxes('Parent',fig,...
    'Position',[15 10 400 300]);

for kk = 1:length(profiles{pp})
    p = plot(ax,profiles{pp}{kk}.pos,profiles{pp}{kk}.z);
    hold(ax,'on')
    pDates{kk} = num2str(profiles{pp}{kk}.date);
end
title(ax,profNums(pp))
ylabel(ax,'Elevation (m NAVD88)','rotation',90)
xlabel(ax, 'Distance (m)')
legend(ax,pDates,'Location','southeast')

dd = uidropdown(fig,...
    'Items',itemNames,...
    'Value', itemNames{1},...
    'Position',[20 310 120 22],...
    'ValueChangedFcn',@(dd,event) selection(dd,p,profiles,ax,profNums));



btn = uibutton(fig,'Position',[150 310 80 22],'Text','export data',...
    'ButtonPushedFcn', @(btn,event) btnGotPushed(profNums,profiles,dd));
end

function selection(dd,p,profiles,ax,profNums)
temp = extractBetween(dd.Value,'e ',' (');
pp = find(profNums==str2double(temp{1}));

cla(ax);
    for kk = 1:length(profiles{pp})
        p = plot(ax,profiles{pp}{kk}.pos,profiles{pp}{kk}.z);
        hold(ax,'on')
        pDates{kk} = num2str(profiles{pp}{kk}.date);
    end
title(ax,profNums(pp))
legend(ax,pDates,'Location','southeast')

setappdata(dd,'p',pp)
end

function btnGotPushed(profNums,profiles,dd)
m=[];   
pp = getappdata(dd,'p');
for kk = 1:length(profiles{pp})
        m = [m;...
            ones(length(profiles{pp}{kk}.z),1)*profNums(pp),... % profile number
            ones(length(profiles{pp}{kk}.z),1)*profiles{pp}{kk}.date,... % profile date
            profiles{pp}{kk}.pos,...                            % profile positon distance
            profiles{pp}{kk}.z];                                % profile elevation
end
t = array2table(m);
t.Properties.VariableNames(1:4) = {'ProfNum','date','distance','elevation'};  
[outfile, outpath] = uiputfile('*.csv',...
    'Choose output location and filename',...
    strcat('Profile_',num2str(profNums(pp),'%03.f'),'_aligned'));
outfile = fullfile(outpath,outfile);
writetable(t,outfile)
end