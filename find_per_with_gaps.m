function out=find_per(inx,iny,D)

%find_per; SWASH processing: creates a reference line for shoreline calculations
%
%This is an mfile to:
%   1. read in a line, 
%   2. find points on the line at a regular spacing, and 
%   3. calculate a perpendicular line at each regularly spaced point.
%       (known as a 'reference line perpendicular').
%I keep track of these reference line perpendiculars
%  through their slope and y intercept (as in the equation y=mx+b. Using
%  the x,y of the reference point, we know all the variables in this eqn)
%
%format:
%
% out=find_per(inx,iny,spacing);
%
%  inx, iny  are vectors of the input line x,y values
%  spacing   is the desired distance between the output lines.
%
% out is a matrix of:
%       1        2       3 4  5
%     [slope y_intercept x y tag]
%      
%      slope is the slope of the reference line perependicular 
%      y_intercept is the y intercept of the reference line perpendicular 
%      x,y is the location of the reference point, thru wich the reference line
%           perpendicular goes
%      tag is the number of the reference point
%
%   caution:   the program does not work well if the distance between the
%               given points is on the order of the desired spacing
%
%afarris@usgs.gov
%afarris@usgs.gov 2017march06 made changes so that large gaps are not
%                 filled with points

% this prgm was started by amy farris sept 26, 96 for nantucket gps work
% altered 8/1/97 for north carolina, the endpoints for the origenal line
%  used are:[421938,4038531 460899,3943473]
%
% First find x,y points along data line at designated spacing:
%  done by finding the distance between two points on the input line, and 
%  finding the angle their connecting line makes with the x-axis: 
%  tan(alpha)=opp/adj.  Then, using the required distance (=D), the x,y of 
%  the first point (x1,y1) and the angle, trigonomic eqs
%  (which equations are used depends on the relative postions of the two pts)
%  are used to solve for x,y of the point D away from x1,y1 (sample eqn:
%   cos(alpha)=(x2-x1)/hyp is solved for x2 (where hyp=required dist)).
%  Then, depending on whether the distance to the next point is > or < the
%  required distance, it moves on to the next point to repeat the procedure.
%  (this is basically the program equal.m)
%
% The perpendicular line is found by calculating the slope of the two lines
%  defined by (x(i-1),y(i-1),x(i),y(i)) and (x(i),y(i),x(i+1),y(i+1)).  
%  These slopes are averaged
%  together and the line perpendicular to point (x(i),y(i)) has
%  a slope that is the negative inverse of the averaged slope.

% tag is a vector the same length as the equally spaced points;  it numbers
%  each point so that it will be easy to keep track of the points.

% firstPoint and lastPoint keep track of big gaps in reference line
% they are the 'first point' and 'last point' in each segment
% where each segment is seperated from other segmentst by large gaps

firstPoint = [];
lastPoint = [];
x1=inx(1); y1=iny(1);
i=2;
done=0;
newx=x1; newy=y1;    
a = 1;
while done == 0 
      if x1<inx(i) & y1<iny(i)
         alpha = atan((iny(i)-y1)/(inx(i)-x1));
         dist2pt=(inx(i)-x1)/cos(alpha);
         xd=x1+D.*cos(alpha);
         yd=y1+D.*sin(alpha);
      elseif x1<inx(i) & y1>iny(i)
         alpha = atan((inx(i)-x1)/(y1-iny(i)));
         dist2pt=(inx(i)-x1)/sin(alpha);
         xd=x1+D.*sin(alpha);
         yd=y1-D.*cos(alpha);
      elseif x1>inx(i) & y1<iny(i)
         alpha = atan((iny(i)-y1)/(x1-inx(i)));
         dist2pt=(x1-inx(i))/cos(alpha);
         xd=x1-D.*cos(alpha);
         yd=y1+D.*sin(alpha);
      elseif x1>inx(i) & y1>iny(i)
         alpha = atan((x1-inx(i))/(y1-iny(i)));
         dist2pt=(x1-inx(i))/sin(alpha);
         xd=x1-D.*sin(alpha);
         yd=y1-D.*cos(alpha);
      elseif x1==inx(i) & y1<iny(i)
         dist2pt=iny(i)-y1;
         xd=x1;
         yd=y1+D;
      elseif x1==inx(i) & y1>iny(i)
         dist2pt=y1-iny(i);
         xd=x1;
         yd=y1-D;
      elseif x1<inx(i) & y1==iny(i)
         dist2pt=inx(i)-x1;
         xd=x1+D;
         yd=y1;
      elseif x1>inx(i) & y1==iny(i)
         dist2pt=x1-inx(i);
         xd=x1-D;
         yd=y1;
      end         
      newx = [newx; xd];
      newy =[newy; yd];
      x1=xd;  y1=yd;
      a = a+1;
%      if dist2pt - D < 0.3.*D  
%       if dist2pt - D < D  
%          i=i+1;         % move on to next data pt
%       end
        doneI = 0;
      if i >= length(inx)
         done=1;
      end
      % I added this bit, b/c results were weird if desited profile spacing
      % was similar or greater than input data spacing.
      % afarris@usgs.gov 2017March02
        while doneI == 0
            dd = dist(xd,yd,inx(i),iny(i));
            if dd >= D
                % next data point is further than desired spacing, 
                % OK to continue normally, keep 'i' the same
                doneI = 1;
            else
                % next data point is too close, move on to next one
                i = i+1;
            end
            if i >= length(inx)
                 doneI = 1;
            end
        end
        % also check if there is too big of gap, 
        if  i < length(inx)
            dd2 = dist(xd,yd,inx(i),iny(i));
            if dd2 > 100 
                % there is a big gap; skip it and restart line on other side
                xd = inx(i);     yd = iny(i);
                x1 = inx(i);     y1 = iny(i);
                newx = [newx; xd];
                newy =[newy; yd];
                i = i+1;
                % keep track of gaps
                lastPoint = [lastPoint; a];
                firstPoint = [firstPoint; a+1];
                a=a+1;
            end
        end
end

% now that we have evenly spaced x,y; time to find perpendicular lines at
%   the points newx,newy
out=zeros(length(newx)-2,4);
a=0;
for i=2:(length(newx)-1)
   a=a+1;
   % fisrt find slope of line on either side of point   
   ma = (newy(i-1)-newy(i))/(newx(i-1)-newx(i));
   mb = (newy(i)-newy(i+1))/(newx(i)-newx(i+1));
   % now calc average slope
   ave_m = (ma+mb)/2;
   % now calc slope of perpendicular line:
   new_m = -1./ave_m;
   % now calc y-intersecpt of perpendicular line:
   b = newy(i) - (new_m).*newx(i);
%   out(a,:)=[new_m b newx(i) newy(i) tag(i)];
   out(a,:)=[new_m b newx(i) newy(i)];
end
% I added the following bit so that first and last points were not skipped
% afarris@usgs.gov 2017March02
% add first and last point, use slope from adjacent profile
first_m = out(1,1);
first_b = newy(1) - (first_m).*newx(1);
last_m = out(end,1);
last_b = newy(end) - (last_m).*newx(end);
%tag = 2:(length(out)+1);
tag = 2:(size(out,1)+1);

% put everything together
out = [first_m first_b newx(1) newy(1) 1;
       out tag';
       last_m last_b newx(end) newy(end) tag(end)+1];

% The slope for the first and last point in a segment are wrong. 
% (since I am averaging in a slope from the other side of a gap)
% Here I'll just use the slope of the adjecent point

for i = 1:length(firstPoint)
    % get slope from next point
    tempSlope = out(firstPoint(i)+1,1);
    % calc correct y intercept using it ( b = y -mx)
    temp_b = out(firstPoint(i),4) - (tempSlope).*out(firstPoint(i),3);
    % corect data
    out(firstPoint(i),1) = tempSlope;
    out(firstPoint(i),2) = temp_b;
end

for i = 1:length(lastPoint)
    % get slope from next point
    tempSlope = out(lastPoint(i)-1,1);
    % calc correct y intercept using it ( b = y -mx)
    temp_b = out(lastPoint(i),4) - (tempSlope).*out(lastPoint(i),3);
    % corect data
    out(lastPoint(i),1) = tempSlope;
    out(lastPoint(i),2) = temp_b;
end










