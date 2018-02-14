% CREATE_FILE_FOR_DSAS_MARSH makes a marsh shoreline data file to be read into DSAS
%
% This code does 4 things:
%  1. loads all the shoreline files (saved by the marsh shoreline code)
%  2. loads the 'bad data' files and uses them to:
%     a.  remove bad data
%     b.  get rid of NaNs in gaps that should be linearly interpoated over
%  4.  Writes out the result as several formats: kml, simple text, old
%       style text, shapefile
%
% Before this code is run, you need to have run marsh_clenup_gui.  It will
% produce a file with the date in the filename, surrounded by '_' in the 
% format of YYYYMMDD.
%
% In ordr for this code to work, you also need to have made a file called
% fileNames.txt.  This file continas the names of all the data files in 
% the region.  It has two columns, the first is the name of the data file,
% the second is the name of the cleanup file produced by marsh_cleanup_gui.
%
% this code is set up to use the smoothed version of the shoreline, but
% you can change this by changing the commenting of the line defining 'temp'
%
% Arc is not very happy with the shapefiles prduced by Matlab.  It is best
% to load the shapefile into Global Mapper and then output it again.
%
% The code is set up to output the data in UTM so it can be compared to the
% other shoreliens, but you can change this by commenting and un-commenting
% parts of the end of the code.
%

% afarris@usgs.gov 2018feb14  now ouputs data as shapefile
% afarris@usgs.gov 2017sep12  made version for marsh shoreline
% afarris@usgs.gov 2014mar10  now can give correct date to each profile
% afarris@usgs.gov 2013may03  added ability to pick bias file
% afarris@usgs.gov 2008june  add date for each profile
% afarris@usgs.gov 2008may ...general updates 
% afarris 2007oct 
% afarris 2007june18 preliminary test version "data4emily.m"

origDir = pwd;

stateDirName = uigetdir(pwd,'Select directory containing shoreline data');

state = input('Enter 2 letter state abbreviation: ','s' );
 
region = input('Enter 2 letter region designation (use "na" if none:  ','s');
region = lower(region);

cd(stateDirName)

%% first get names of files to load
 
% For simplicity's sake just demand that there is a file containing the
% names of the files in the order they need to be read in.  Also demand
% that this file be called 'fileNames.txt'  It has 2nd column which is 
% the name of the cleanup files 
 
% This is set up to load all the files for a given region at one time.

fname = 'fileNames.txt';

% load file that contians the names of data files to load
fid =fopen(fname,'r');
if fid < 3
    error('could not find file')
end
C = textscan(fid,'%s %s');
fclose(fid)
% C{1} is the filenames of shoreline data
% C{2} is the filename of the cleanup file

%% now load data

% preallocate matricies to hold all the data
data=[];
dateArray=[];
allTrNum=[];
lastTrNum = 0;

numFiles = size(C{1},1);
for i = 1: numFiles% first load shoreline file
    % the following line loads the file.  in it are several variables
    % we need marshSmoothXYall
    cellfun(@load,C{1}(i));
    disp(C{1}(i))
    temp = marshSmoothXYall;
%    temp = marshXYZmAll(:,1:2);
    
    % now load cleanup file
    % loads variables prfiles2delet and fillGap
    cellfun(@load,C{2}(i));
    
    % get date from filename
    nameTemp = char(C{2}(i));
    ff = strfind(nameTemp,'_');
    junk = nameTemp(ff(1)+1:ff(2)-1);
    fancyDateStr = [junk(5:6),'/',junk(7:8),'/',junk(1:4)];
    foo2=repmat(fancyDateStr,size(temp,1),1);
    tempDate = cat(1,dateArray,foo2);
    
    figure
    plot(temp(:,1),temp(:,2),'r*')
    hold on
    % use profiles2dete to delete priofiles
    for ii=1:length(profiles2delete)
        foo = str2num(cell2mat(profiles2delete(ii)));
        f = find(ptNum == foo);
        if ~isempty(f)
            temp(f,1:2) = [NaN NaN ];
        end
    end
    plot(temp(:,1),temp(:,2),'k.-')
    title(nameTemp)
    
    % now add data to matrix that holds all the data
    data = [data;  NaN NaN; temp];
    % ptNum is missing the last point
    ptNum = [ptNum; ptNum(end)+1];
    % add to vecotr of all transect numbers
    allTrNum = [allTrNum; ptNum + lastTrNum];
    lastTrNum = allTrNum(end);
    dateArray = [dateArray; fancyDateStr ; tempDate];
end 

%% make output data files

%% write out data
zone = input('what is the UTM zone?  ');

% create output filename 
if strcmp(region,'na') == 1
    guess = {[state,'_marsh_shoreline_for_DSAS.txt']};
else
    guess = {[state,region,'_marsh_shoreline_for_DSAS.txt']};
end

% check to make sure that this filename is OK
prompt = {'Enter an output file name; '};
dlg_title = 'Output file name';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines,guess);
outName = answer{1};

fid=fopen(outName,'wt');
% This file needs a header that gives the UTM zone
fprintf(fid,'%% This file contains lidar marsh shoreline data.  The columns are: \n');
fprintf(fid,'%%  x  y  date \n');
fprintf(fid,'%%  x,y are UTM (m) easting, northing; zone %2.0f',zone);
fprintf(fid,'   \n');
fprintf(fid,'%%  \n');
for i = 1:length(data)
    fprintf(fid,'%9.3f %10.3f %10s \n',data(i,:)',dateArray(i,:));
end
fclose(fid)

%% create a file I can load
% create output filename 

% create output filename 
if strcmp(region,'na') == 1
    name2 = [state,'_marsh_shoreline_for_DSAS_easy.csv'];
else
    name2 = [state,region,'_marsh_shoreline_for_DSAS_easy.csv'];
end

out = [data(~isnan(data(:,1)),:) allTrNum(~isnan(data(:,1)))];
outDate = dateArray(~isnan(data(:,1)),:);

fid=fopen(name2,'wt');
fprintf(fid,'x, y, trNum, Date_  ');
fprintf(fid,'   \n');
for i = 1:length(out)
    fprintf(fid,' %9.3f , %10.3f, %9.3f, %10s\n',out(i,:)',outDate(i,:));
end
fclose(fid)

%% create kml
if strcmp(region,'na') == 1
    name3 = [state,'_marsh_shoreline_cleaned'];
else
    name3 = [state,region,'_marsh_shoreline_cleaned'];
end

[lonMarsh, latMarsh] = utm2ll(out(:,1),out(:,2),19);
kml_line(lonMarsh, latMarsh,name3)

cd(origDir)

%% create shapefile




% UTM
[Xcells,Ycells] = polysplit(data(:,1),data(:,2));
n = length(Ycells)
for i = 1: n
    [TracksSplit(i).X] = Xcells{i};
    [TracksSplit(i).Y] = Ycells{i};
    bb = [min(Xcells{i}) min(Ycells{i}); max(Xcells{i}) max(Ycells{i})];
    TracksSplit(i).BoundingBox =  bb;
end
name4 = [name3,'UTM'];

% geographic
% [lon, lat] = utm2ll(data(:,1),data(:,2),19);
% [latcells,loncells] = polysplit(lat,lon);
% n = length(loncells)
% for i = 1: n
%     [TracksSplit(i).Lon] = loncells{i};
%     [TracksSplit(i).Lat] = latcells{i};
%     bb = [min(loncells{i}) min(latcells{i}); max(loncells{i}) max(latcells{i})];
%     TracksSplit(i).BoundingBox =  bb;
% end
% name4 = [name3,'GEO'];

[TracksSplit(1:n).Geometry] = deal('Line');
[TracksSplit(1:n).Name] = deal('Marsh edge');
[TracksSplit(1:n).Date_] = deal(dateArray(1,:));
% if dateArray(1,:) ~= dateArray(end,:)
%     % it looks like there is more than one date.
%     keyboard
% end
msgbox('The date listed for every point will be the date of the first point')



shapewrite(TracksSplit,name4)


