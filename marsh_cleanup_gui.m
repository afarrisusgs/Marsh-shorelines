function shoreline_cleanup_gui

% SHORELINE_CLEANUP_GUI interactivly remove bad shorelines and allow gap filling
%
% This code is intented to be run after the shoreline code is run and the
% normal handchecking is done.  This code needs to be run before
% bias_calc.m  It has two main purposes:  1. to delete shoreline points that
% look wrong and 2. to mark which gaps in shoreline position can be filled .
% in DSAS. This code loads the SL_*.mat file from the shoreline code and the
% profile data from which the shorelines were calculated.  It saves a mat
% file with two variables:  profiles2delete and fillGap.
%  Profiles2delete is a vector of profile numbers that need to be deleted
% fillGap is a matrix with 2 columns.  The first is all the profile numbers
% in the SL data file, the second is a flag that is either 0 or 1.  O means
% don't fill gap, 1 means fill gap.
%
% After this code is run, you can run "shoreline_check.m" to check the
% results (it loads the mat file and plots them).
%
% the profiles in 'profiles2delete' are not necessarily in order

%afarris@usgs.gov 2007
%afarris@usgs.gov 2009oct16 finally generalized y-tick labels on colobar


%%  load data files and set up variables

% this code started out very simple.  The problem is that the user might 
% need to see some profiles from the adjacent tiles in order to make 
% descisions about gaps at the edges of tiles. This unfortunately 
% complicates the code, and I couldn't find a way around it.

% Start by loading shorelines and profile data for the main tile of
% interest:
% note that the shoreline data file is a mat file and when it is loaded a
% variable named SL_data is placed in the workspace.  The profile data is
% an ascii file and therefore when it is loaded it can be assigned a 
% variable name.
%
addpath('C:\mfiles\others\riches\')
addpath('C:\mfiles\plotting\')

%% get marsh shoreline
[marshFileName, marshPathName] = uigetfile('*.mat','Click on marsh data file.');
disp(['you have chosen ', marshFileName])
% load marsh shoreline:
load([marshPathName, marshFileName]);                
mainName = marshFileName;

ddate = inputdlg('Date of survey (format YYYYMMDD): ','Date?',1)

% create output name
f=findstr(marshFileName,'.');
outName = [marshPathName,marshFileName(1:f-1),'_',cell2mat(ddate), '_cleanup.mat'];
handles.outName = outName;
%TODO check if I'll be overwriting a file


%% load las data

origDir = pwd;
dirName = uigetdir(pwd,'Click on the directory containing the las data');
cd(dirName)
disp(['looking in ', dirName])

dd = dir;

% load data
xAll = [ ];
yAll = [ ];
zAll = [ ];

for i = 1:length(dd)
    if ~dd(i).isdir
        [pathstr,name,ext] = fileparts(dd(i).name) ;
        if isequal(ext,'.las')
            disp(dd(i).name)
            A = LASreadAll(dd(i).name);
            xAll = [xAll; A.x ];
            yAll = [yAll; A.y ];
            zAll = [zAll; A.z ];

        end
    end
end

clear A
cd(origDir)

%%
% This code was written for shorelines.  The height of a shoreline is at
% MHW.

MHW = nanmean(mhwXYZmAll(:,3));
MTL = nanmean(mtlXYZmAll(:,3));
if isnan(MHW)
    MHW = input('Please enter the approximate MHW:  ')
end
if isnan(MTL)
    MTL = input('Please enter the approximate MTL :  ')
end


% work on colormap, this is based on MHW shoreline work, but since the
% shoreline is not determined by and eleveatoion, I don't force a color
% jump at MHW, I just set upper and lower bounds
zmin = MTL - 0.4;
zmax = MHW + 0.3;

% in order to make everything faster, don't plot all of the lidar data
f = find(zAll >  MTL - 0.4 & zAll<MHW + 0.3);
profileX = xAll(f);
profileY = yAll(f);
profileZ = zAll(f);

shoreNums = ptNum;

shoreX = marshSmoothXYall(:,1);
shoreY = marshSmoothXYall(:,2);


% set up vectors
prNumFirst = shoreNums(1);
% FIXME what if shoreNums(1) is a NaN?
prNumLast = shoreNums(end);
if isnan(prNumLast) 
    % last point is a NaN, find last point that isn't a NaN
     f=isnan(shoreNums);
     ii=find(f==0,1,'last');
     prNumLast = shoreNums(ii);
     %check if it worked:
     if isnan(prNumLast) 
         % it didn't work
         error('trouble with profile numbers (lin 200)')
     end
end
     
allPrNums = prNumFirst:prNumLast;

f=find(isnan(shoreX)==0);
shoreXnoNans = shoreX(f);
shoreYnoNans = shoreY(f);
shoreNumsNoNans = shoreNums(f);

%interpolate
fakeX = interp_rich(shoreNumsNoNans,shoreXnoNans,allPrNums);
fakeY = interp_rich(shoreNumsNoNans,shoreYnoNans,allPrNums);
prNumsAddNan = allPrNums;

% I want to create a vector of 0s and 1s.  0 means don't fill a gap
% and 1 means do fill gap
flag = zeros(size(allPrNums,2),1);

% load colormap 
% load colormapMHWjump20


% Put importatnt data in the 'handles' structure so it can be acessed by
% all the callbacks
handles.shoreX = shoreX;
handles.shoreY = shoreY;
handles.shoreNums = shoreNums;

handles.fakeX = fakeX;
handles.fakeY = fakeY;
handles.allPrNums = allPrNums;

handles.flag = flag;
% handles.padded = padded;

%save data to make plots:
handles.MHW = MHW;
%profile data
%handles.DAP = profiles(:,3);
handles.z = zAll;
% handles.profilePrNum = profiles(:,6);
%shoreline dataones(size(handles.shX));
handles.shX = shoreX;
handles.zdatum = MHW;
handles.shorePrNum = ptNum;

% the following line was changed on 2007sep13
% Now after this code has been run and its ouput file is loaded, you can 
% test to see if any profiles were deleted by typing:
%   isempty(profiles2delete)
% This method doesn't work with ouput from previous version of the code
handles.profiles2delete = [];


%% lay out gui
handles.figure=figure('position',[100 80 620 820]);

set(handles.figure,'toolbar','figure')

% Now add all the buttons; they are organized in columns:

%-------------------------------------------------
% First column -- find profile number
% add button for user to press to get profile number
% (of most recently clicked point)
handles.buttonGetPrNum = uicontrol('style','pushbutton','string',...
    'Find profile number',...
    'position',[30 150 100 30],'parent',handles.figure,'fontSize',8);
% add box to disply profile number
handles.textBoxShowPrNum = uicontrol('Style','text','string',...
    'Profile #','position',[20 100 120 40]);
% add a button to make plot
handles.buttonPlotPr = uicontrol('style','pushbutton','string',...
    'Plot profile ',...
    'position',[30 50 100 30],'parent',handles.figure,'fontSize',8);

%-------------------------------------------------
% Second column -- delete shorelines
% add button for user to press to delete shoreline
handles.buttonDeleteSh = uicontrol('style','pushbutton','string',...
    'Delete Shoreline',...
    'position',[180 150 100 30],'parent',handles.figure,'fontSize',8);
% add text to explain list box of shorelines to delete
handles.textBoxShowPr = uicontrol('Style','text','string',...
    'Profiles to delete:','position',[180 130 100 15]);
% add list box that shows the shorelines to delete
handles.listBoxDeletedSh = uicontrol('style','listbox',...
    'position',[180 40 100 80]);
% add button to allow user to remove the highlighted element in the listbox
handles.buttonRemoveSh = uicontrol('style','pushbutton',...
    'position',[180 5 100 20],'string','Remove Selection',...
    'backgroundColor','y');

%-------------------------------------------------
% Third column -- fill gaps
% add button for user to press to allow filling of gap in shoreline
handles.buttonFillGap = uicontrol('style','pushbutton','string',...
    'Fill Gap',...
    'position',[320 150 80 30],'parent',handles.figure,'fontSize',8);
% add text to explain list box of gaps to fill
handles.textBoxShowGap = uicontrol('Style','text','string',...
    'Gaps to fill:','position',[320 130 100 15]);
% add list box that shows the gaps to fill
handles.listBoxFilledGaps = uicontrol('style','listbox',...
    'position',[320 40 100 80]);
% add button to allow user to remove the highlighted element in the listbox
handles.buttonRemoveGap = uicontrol('style','pushbutton',...
    'position',[320 5 100 20],'string','Remove Selection',...
    'backgroundColor','y');

%-------------------------------------------------
% Fouth column -- quit
% add button for user to press to quit w/o saving
handles.buttonQuit = uicontrol('style','pushbutton','string',...
    'Quit',...
    'position',[450 150 80 30],'parent',handles.figure,'fontSize',8);
% add button for user to press to quit and save
handles.buttonQuitSave = uicontrol('style','pushbutton','string',...
    'Quit and Save',...
    'position',[450 100 80 30],'parent',handles.figure,'fontSize',8);

% add text telling name of file being displayed
handles.textBoxShowPr = uicontrol('Style','text','string',...
    'File being displayed:','position',[450 25 150 18]);
handles.textBoxShowPr = uicontrol('Style','text','string',...
    mainName,'position',[430 10 180 18]);

%% Now add plot to gui
handles.axes = axes('Parent',handles.figure,'Position',[.12 .3 .77 .65],...
     'units','pixels','box','on');

% handles.axes = axes('Parent',handles.figure,'Position',[.12 .35 .77 .60],...
%cdots_amy(xAll,yAll,xAll,zmin,zmax,cmapMHWjump20,handles.axes)
cdots_amy(profileX,profileY,profileZ,zmin,zmax,jet,handles.axes)
hold on
% keep axis limits in case of really wild shoreline points:
axis(axis)
nanplot(shoreX,shoreY,'color',[.75 .75 .75])
plot(shoreX,shoreY,'k-','linewidth',2)
p2 = plot(fakeX,fakeY,'o','linewidth',2,'markerEdgeColor',...
    'k','markerFaceColor','k','markerSize',6);
p1 = plot(shoreX,shoreY,'o','linewidth',2,'markerEdgeColor',...
    'k','markerFaceColor','w','markerSize',6);
axis equal
title('Colored dots: Profile data , white circles: shoreline location')
xlabel('UTM Easting (m)')
ylabel('UTM Northing (m)')


% make colorbar
% calc the change in height per bin in colormap (there are 64 bins)
zinc=(zmax-zmin)/64;
% the ticks will be at
yt = [10 20 30 40 50 60];
% calculate the labels for these ticks
% b tells us where we are in the cell 'ytl_string'
b=1;
for a= yt
    ytl(b) = zmin + a*zinc; 
    ytl_string(b) = {num2str(ytl(b),'%3.1f')};
    b=b+1;
end
ytl_string{2}= 'MHW';

cb = colorbar('ytick',yt,'yticklabel',ytl_string);
title(cb,'Height (m)')
colormap(jet);

%% define callbacks
set(handles.buttonGetPrNum, 'callback',{@buttonGetPrNum_callback, handles});
set(handles.buttonPlotPr, 'callback',{@buttonPlotPr_callback, handles});
set(handles.buttonDeleteSh, 'callback',{@buttonDeleteSh_callback, handles});
set(handles.buttonRemoveSh, 'callback',{@buttonRemoveSh_callback, handles});
set(handles.buttonFillGap, 'callback',{@buttonFillGap_callback, handles});
set(handles.buttonRemoveGap, 'callback',{@buttonRemoveGap_callback, handles});
set(handles.buttonQuit, 'callback',{@buttonQuit_callback, handles});
set(handles.buttonQuitSave, 'callback',{@buttonQuitSave_callback, handles});

guidata(handles.figure,handles)

test(1) = cellstr('In this plot the lidar profile data is plotted as colored dots.');
test(2) = cellstr('The shoreline data is plotted on the top of the profile data as white dots');
test(3) = cellstr('For some profiles, no shoreline position was calclated. ');
test(4) = cellstr('These gaps were linearly interpolated over and shown as black dots.  ');
test(5) = cellstr('If you want DSAS to fill these gaps you need to click on each black dot,');
test(6) = cellstr('then on the "Fill Gap" button.  The point you selected will turn grey. ');
test(7) = cellstr('If you want to delete a shoreline point. Click on it, then click on "Delete profile".');
test(8) = cellstr('The point you selected will have a red x on it.  ');
test(9) = cellstr('You can zoom and pan this figure, but you must exit these features before you click on a profile.');
test(10) = cellstr('Note that all white button work on the last profile the user cliked on.');
test(11) = cellstr('While the yellow buttons work on the highlighted element in their ');
test(12) = cellstr('respective listbox.');

h = msgbox(test);

end % shoreline_cleanup_gui


%% callbacks
% ------------callback for button to find profile number-------------
function buttonGetPrNum_callback(hObject,event_data,handles)

handles = guidata(hObject);

% get the profile number from the user's last click on the plot 
p = get(handles.axes,'currentPoint');
d=sqrt((handles.fakeX - p(1,1)).^2 + (handles.fakeY-p(2,2)).^2);
% I use the 'fake' values so that the user can get the profile numbers of
% profiles with no shoreline data
[i,j] = min(d);
chosenProfileNum = handles.allPrNums(j);
set(handles.textBoxShowPrNum,'String',['Profile # is ',num2str(chosenProfileNum)]);

% TODO label this shoreline on map view plot?

end % buttonGetPrNum_callback

% ------------callback for button to plot profile -------------
function buttonPlotPr_callback(hObject,event_data,handles)

handles = guidata(hObject);

% get the profile number from the user's last click on the plot 
p = get(handles.axes,'currentPoint');
d=sqrt((handles.fakeX - p(1,1)).^2 + (handles.fakeY-p(2,2)).^2);
% I use the 'fake' values so that the user can get the profile numbers of
% profiles with no shoreline data
[i,j] = min(d);
chosenProfileNum = handles.allPrNums(j);
% find rows in profile data with the selected profile number
fPr = find(handles.profilePrNum==chosenProfileNum);
% find row in shoreline data with the selected profile number
fSh = find(handles.shorePrNum==chosenProfileNum);
% now make plot
figure
plot(handles.DAP(fPr),handles.z(fPr),'.')
hold on
% add shoreline position with error
plot(handles.shX(fSh),handles.zdatum(fSh),'r*')
% errmin = SL_data(shRowNum,6)-.5*SL_data(shRowNum,22);
% errmax = SL_data(shRowNum,6)+.5*SL_data(shRowNum,22);
% plot([errmin errmax],[SL_data(shRowNum,22) SL_data(shRowNum,22)],'r')
% add MHW
ax = axis;
axis manual
plot([ax(1) ax(2)],[handles.MHW handles.MHW],'k')
if isnan(handles.shX(fSh))
    % find row number in handes.ref for this profile number
    pp = find(handles.ref(:,5)==chosenProfileNum);
    dd = dist(handles.fakeX(fSh),handles.fakeY(fSh),handles.ref(pp,3),handles.ref(pp,4))
    % FIXME is there a way to get a sign for dd?
    plot([-1*dd -1*dd],[ax(3) ax(4)],'k')
    plot([dd dd],[ax(3) ax(4)],'k')
end
xlabel('Distance from profile origin (m)')
ylabel('Height (m)')
title(['Profile number ',num2str(chosenProfileNum)]);

end % buttonPlotPr_callback


% ----------------Callback for button to delete shoreline -------------------
% Executes after button press
function buttonDeleteSh_callback(hObject, eventdata, handles)
% hObject    handle to edit_profile_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% get the profile # of the user's last click on the plot
p = get(handles.axes,'currentPoint');
d=sqrt((handles.shoreX - p(1,1)).^2 + (handles.shoreY-p(2,2)).^2);
[i,j] = min(d);
chosenProfileNum = handles.shoreNums(j);

% check to see if a profile number is already in list box
junk = get(handles.listBoxDeletedSh,'string');
if size(junk,1) > 0
    % add this profile number to a list
    temp = [cellstr(num2str(chosenProfileNum)); junk];
else
    % make the list starting with this number
    temp = cellstr(num2str(chosenProfileNum));
end
set(handles.listBoxDeletedSh,'string',temp);
handles.profiles2delete = temp;
beep 

%add to plot
f =find(handles.shoreNums==chosenProfileNum);
plot(handles.shoreX(f),handles.shoreY(f),'+','linewidth',2,'markerEdgeColor',...
    'r','markerFaceColor','r','markerSize',6);

guidata(hObject, handles);
end % buttonDeletedSh_callback

% ------------callback for button to remove highlighted selection-------------
function buttonRemoveSh_callback(hObject, eventdata, handles)
% hObject    handle to edit_profile_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% get users selection
selectedIndex = get(handles.listBoxDeletedSh,'value');
allStrings = get(handles.listBoxDeletedSh,'string');
temp = allStrings;
% remove it
temp(selectedIndex)=[];
% re-dispaly 
set(handles.listBoxDeletedSh,'string',temp,'value',1);
handles.profiles2delete = temp;

chosenProfileNum = cell2mat(allStrings(selectedIndex));
f =find(handles.shoreNums==str2num(chosenProfileNum));
plot(handles.shoreX(f),handles.shoreY(f),'o','linewidth',2,'markerEdgeColor',...
    'k','markerFaceColor',[1 1 1],'markerSize',6)

guidata(hObject, handles);
end % buttonRemoveSh_callback

% -------------------callback for button to identify gap -----------------
function buttonFillGap_callback(hObject,event_data,handles)

handles = guidata(hObject);

%get data out of handles
fakeX=handles.fakeX;
fakeY=handles.fakeY;
allPrNums=handles.allPrNums;

%get number of profile that was clicked on
p = get(handles.axes,'currentPoint');
d=sqrt((fakeX - p(1,1)).^2 + (fakeY-p(2,2)).^2);
[i,j] = min(d);
chosenProfileNum = allPrNums(j);
f =find(handles.allPrNums==chosenProfileNum);
%add to plot
p2 = plot(fakeX(f),fakeY(f),'o','linewidth',2,'markerEdgeColor',...
    'k','markerFaceColor',[.5 .5 .5],'markerSize',6);

% check to see if a profile number is already in list box
junk = get(handles.listBoxFilledGaps,'string');
if size(junk,1) > 0
    % add this profile number to a list
    temp = [cellstr(num2str(chosenProfileNum)); junk];
else
    % make the list starting with this number
    temp = cellstr(num2str(chosenProfileNum));
end
set(handles.listBoxFilledGaps,'string',temp);
beep 

%update handles
handles.fakeX = fakeX;
handles.fakeY = fakeY;
handles.flag(f) = 1;

guidata(hObject,handles)
end % buttonFillGap_callback

% ------------callback for button to remove highlighted selection-------------
function buttonRemoveGap_callback(hObject, eventdata, handles)
% hObject    handle to edit_profile_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% get users selection
selectedIndex = get(handles.listBoxFilledGaps,'value');
allStrings = get(handles.listBoxFilledGaps,'string');
temp = allStrings;
chosenProfileNum = cell2mat(allStrings(selectedIndex));
% remove it
temp(selectedIndex)=[];
% re-dispaly 
set(handles.listBoxFilledGaps,'string',temp,'value',1);
%change data
f =find(handles.allPrNums==str2num(chosenProfileNum));
handles.flag(f) = 0;
%add to plot
p2 = plot(handles.fakeX(f),handles.fakeY(f),'o','linewidth',2,'markerEdgeColor',...
    'k','markerFaceColor',[0 0 0],'markerSize',6);

guidata(hObject, handles);
end % buttonRemoveGap_callback

% ---------------------callback for quit button ------------------------
function buttonQuit_callback(hObject,event_data,handles)
close(handles.figure)
end %buttonQuit_callback

% ---------------------callback for quit and save button ------------------
function buttonQuitSave_callback(hObject,event_data,handles)

handles = guidata(hObject);

% get name of output file
name = handles.outName;

%get variables to write out
fillGap = [(handles.allPrNums)' handles.flag];
profiles2delete = handles.profiles2delete;

    
% save data
save(name, 'fillGap', 'profiles2delete')
disp('The gap flag was saved in a file called:')
disp(name)

% all done
msgbox('data has been saved.  you can close the gui at will')

end %buttonQuitSave_callback


