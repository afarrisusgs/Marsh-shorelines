% CREATE_FILE_FOR_DSAS_MARSH makes a marsh shoreline data file to be read into DSAS
%
% This code does 4 things:
%  1. loads all the shoreline files (saved by the marsh shoreline code)
%  2. loads the 'bad data' files and uses them to:
%     a.  remove bad data
%     b.  get rid of NaNs in gaps that should be linearly interpoated over
%  4.  Writes out a simple text file of
%           1  2   3      
%           x  y  date]
%
%  x,y are UTM easting and northing in meters.  
% Any NaNs still in the file designate gaps that should NOT be interpolated
% over.  'bad data' files created by "marsh_cleanup_gui.m", the date is
% in the filename, surrounded by '_' in the format of YYYYMMDD.
%
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

numFiles = size(C{1},1);
for i = 1: numFiles% first load shoreline file
    % the following line loads the file.  in it are several variables
    % we need marshSmoothXYall
    cellfun(@load,C{1}(i));
    disp(C{1}(i))
    temp = marshSmoothXYall;
    
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

out = data(~isnan(data(:,1)),:);
outDate = dateArray(~isnan(data(:,1)),:);

fid=fopen(name2,'wt');
fprintf(fid,'x, y, z ');
fprintf(fid,'   \n');
for i = 1:length(out)
    fprintf(fid,'%9.3f, %10.3f,  %10s  \n',out(i,:)',outDate(i,:));
end
fclose(fid)

cd(origDir)