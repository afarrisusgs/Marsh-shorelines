%% set path on laptop
addpath('C:\Projects\MassCZM\newMethod')
addpath('C:\mfiles')
addpath('C:\mfiles\SVN\profile_cutting')
addpath('C:\mfiles\SWASH_processing')
addpath('C:\mfiles\others\m_contrib\trunk\m_map')
addpath('C:\mfiles\others\m_cmg\trunk\RPSstuff')
addpath('C:\mfiles\others\matlab_central')
addpath('C:\mfiles\others\riches\')
addpath('C:\mfiles\plotting\')

%% load las data
% region = input('Enter 2 letter region designation (use "na" if none:  ','s');

origDir = pwd;
dirName = uigetdir(pwd,'Click on the directory containing the las data');
cd(dirName)
disp(['looking in ', dirName])

% cd C:\Users\afarris\Desktop\MassCZM\newMethod\sandwhich
% cd C:\Projects\MassCZM\newMethod\sandwich
% cd C:\Projects\MassCZM\newMethod\chatham

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

if isempty(xAll)
    error('No las files found!')
end

 
%% load contour
% cd C:\Projects\MassCZM\newMethod\sandwich
% cc = csvread('sand_zero_contour.csv',1,0);

[fileName, pathname] = uigetfile('*.csv','Click on the contour file');
disp(['loading file ' fileName])

cc = csvread(fullfile(pathname, fileName),1,0);


%%

% set up mesh that grid will be based on. 
xx = ceil(min(xAll)):1:floor(max(xAll));
yy = ceil(min(yAll)):1:floor(max(yAll));
[X,Y] = meshgrid(xx,yy);
% grid data
Z = griddata(xAll,yAll,zAll, X, Y);


%% pass data to marsh code
% ask user how they want the codes to run

prompt = {'What value do you want to use for the threshold?, Enter 0 for none',...
    'How far offshore of MHW do you want to start looking for the marsch edge?',...
    'What is MHW?', 'What is MTL?'};
dlg_title = 'Code customization';
num_lines = 1;
default_ans = {'3','5','1.22','0'};
answer = inputdlg(prompt,dlg_title,num_lines,default_ans);
threshold = str2num(answer{1});
offset = str2num(answer{2});
MHW = str2num(answer{3});
MTL = str2num(answer{4});
name = input('Name for output file:  ','s');
marshName = [name,'_marsh_t',num2str(threshold),'_o',num2str(offset)]

seg = cc(:,4);
% northern part of sandwich:
% seg = cc(1:2663,4);
% southern part of chatham:
% seg = cc(2203:2980,4);

segU = unique(seg);
mhwXYZmAll = [ ];
marshXYZmAll = [ ];
mtlXYZmAll = [ ];
marshSmoothXYall = [ ];
endXYZall = [ ];
% this method does not keep profile number, but I need to number the points
ptNum = [ ];

for i = 1:length(segU)
% for i = 31 : 57
    f = find(cc(:,4)==segU(i));
    c = cc(f,:);
    % first see how long the segement is
    dd = dist(c(1,5),c(1,6),c(end,5),c(end,6));
    if dd > 20
        % caluclate marsh edge
        [mhwXYZm,mtlXYZm,marshXYZm,marshSmoothXY,endXYZ] = ...
            calc_marsh_edge(X,Y,Z,c,threshold,offset,MHW,MTL);
        % put NaNs between sections or kml looks bad
        gap4 = NaN .* ones(1,4);
        gap3 = NaN .* ones(1,3);
        gap2 = NaN .* ones(1,2);
        % save results from this section in final matrix
        mhwXYZmAll = [mhwXYZmAll; gap4; mhwXYZm];
        marshXYZmAll = [marshXYZmAll; gap4; marshXYZm];
        mtlXYZmAll = [mtlXYZmAll; gap4; mtlXYZm];
        marshSmoothXYall = [marshSmoothXYall; gap2; marshSmoothXY];
        endXYZall = [endXYZall; gap3; endXYZ];
        % number the points (and an extra number for the gap between this
        % section and the next)
        foo =  1 : size(endXYZ,1)+1;
        if isempty(ptNum)
            ptNum = foo';
        else
            ptNum = [ptNum; foo' + ptNum(end)];
        end
    end
end

%% make kml file and save data in mat file.

[lonMarsh, latMarsh] = utm2ll(marshXYZmAll(:,1),marshXYZmAll(:,2),19);
kml_line(lonMarsh, latMarsh,marshName)

[lonBeach, latBeach] = utm2ll(mhwXYZmAll(:,1),mhwXYZmAll(:,2),19);
kml_line(lonBeach, latBeach,[name, '_beach'])

[lonEnd, latEnd] = utm2ll(endXYZall(:,1),endXYZall(:,2),19);
kml_line(lonEnd, latEnd,[name,'_end'])

f= find(marshSmoothXYall(:,1)> 1e5);
[lonS, latS] = utm2ll(marshSmoothXYall(f,1),marshSmoothXYall(f,2),19);
kml_line(lonS, latS,[marshName,'_smooth'])

save(name,'mhwXYZmAll','mtlXYZmAll','marshXYZmAll','marshSmoothXYall','endXYZall','ptNum','MHW')

