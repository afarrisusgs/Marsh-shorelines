function distance = dist(x1,y1,x2,y2)
% DIST calculate the distance between two points or many points
%
%  Usage:
%      distance=dist(x1,y1,x2,y2)
%
%         output:
%            distance:  the distance between two points, or a series of 
%                       points.  It will be the same size as x2 and y2.
%         input:
%            x1,y1:     starting point, a single point or vectors of pts
%            x2,y2:     ending point or vector or matrix of points. 
%                       If vectors, the first point of x2,y2 would usually 
%                       be x1,y1, so that the first distance calculated is zero.
%                       If (x2,y2) are matricies, the distance is found 
%                       columnwise between each value in (x1,y1) and all
%                       the poitns in each column of (x2,y2).
%                       If (x1,y1) are vectors, their length must equal
%                       the number of columns in (x2,y2).
%
%  Examples:
%          % distance between two points
%          distance=dist(x1,y1,x2,y2)
%          % distance between many points and the first point
%          distance=dist(x(1),y(1),x,y)
%           % now show that x2,y2 can be matricies:
%           x1=[0 2 4];
%           y1=[0 2 4];
%           x2=[3 5 7 ;5 7 9 ;7 9 11 ;8 10 12 ];
%           y2=[4 6 8 ;12 14 16 ;24 26 28 ;15 17 19 ];
%           distance=dist(x1,y1,x2,y2)

% afarris@usgs.gov
% jlist 2/17/05:  comments upgraded, and variable 'dist' changed to distance
% afarris 2007july25: changed so that x2,y2 can be matricies

%% check that all input variables are the right sizes
[x1r,x1c]=size(x1);
[y1r,y1c]=size(y1);
if x1r > 1 && x1c > 1
    error('first input argument must be a point or a vector, not a matrix')
end
if y1r > 1 && y1c > 1
    error('first input argument must be a point or a vector, not a matrix')
end

if x1r ~= y1r || x1c ~= y1c
    error('the first 2 input variables must be the same size')
end

[x2r,x2c]=size(x2);
[junk1,junk2]=size(y2);
if junk1 ~= x2r || junk2 ~= x2c
    error('third and fourth variables must be the same size')
end

% if max(x1r,x1c) >1 && max(x1r,x1c) ~= x2c
%     error('length of first 2 input arguments must be the same as the number of columns in the last 2 input arguments')
% end


%% do calculation

% pre-allocate output matrix
distance = NaN * ones(x2r,x2c);

if x1r == 1 && x1c ==1
    % first input arguments are points, not vector
    for i = 1:x2c
        distance(:,i)=sqrt((x1-x2(:,i)).^2 +(y1-y2(:,i)).^2);
    end
else
    % first input arguments are vectors
    for i = 1:x2c
        distance(:,i)=sqrt((x1(i)-x2(:,i)).^2 +(y1(i)-y2(:,i)).^2);
    end
end

