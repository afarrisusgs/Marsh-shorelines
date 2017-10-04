function [mhwXYZm,mtlXYZm,marshXYZm,marshSmoothXY,endXYZ]  = calc_marsh_edge(X,Y,Z,c,threshold,offset,MHW,MTL)

% this code calculates a marsh edge from lidar data
% X,Y,Z are the gridded elevation data
% c is the MTL contour that will be used as a baseline
% name is the name of the region....to be used to name the GE files 

% afarris@usgs.gov 2017may25
% afarris@usgs.gov 2017march01


%% find ref line perpendiculars aka transects
%  I am treading these transects as ephemeral, that is I am not saving them

inx = c(:,5);
iny = c(:,6);
    
perp=find_per_with_gaps(inx,iny,5);

% figure
% plot(inx,iny,'k-*')
% hold on
% plot(perp(:,3),perp(:,4),'r.-')
% axis equal

%% now step thru each transect to find marsh edge

% set up matrix to hold output
% keep (x,y,z) and slope at MHW and at marsh (if there is one)
mhwXYZm = NaN.*ones(size(perp,1),4);
mtlXYZm = NaN.*ones(size(perp,1),4);
marshXYZm = NaN.*ones(size(perp,1),4);
marshSmoothXY = NaN.*ones(size(perp,1),2);
endXYZ = NaN.*ones(size(perp,1),3);

for i = 1 : length(perp)
    % in order to get DEM along these transects, I need the transects in the
    % format of (x,y) pairs. Currently, I have the slope and y intercept of the
    % line.  So I guess I can come up with a series of x values and calculate
    % the y values.  How do I come up with the x values??  Easiest to say xmin
    % = xorg -10, xmax = xorg+10; But if line is very N-S this technique 
    % may not work. So I'll check slope, if line runs N-S, I'll set y
    % values.  You can look at Matlab's docutmentaion for atan2 to see the 
    % meaning of slope (perp(:,1))
    if perp(i,1) > pi || perp(i,1) < -pi
        % line is N-S
%        yt = perp(i,4)-50: 0.2 : perp(i,4)+50;
         yt = perp(i,4)-30: 0.2 : perp(i,4)+30;
        xt = (yt - perp(i,2)) / perp(i,1);
    else
        % line is East-West
%        xt = perp(i,3)-50: 0.2 : perp(i,3)+50;
        xt = perp(i,3)-30: 0.2 : perp(i,3)+30;
        yt = perp(i,1) * xt + perp(i,2);
    end
    % xt, yt are (x,y) points along the transect
    
    % now get elevation from DEM at these locations
    zt = interp2(X,Y,Z,xt,yt);

    % OK, so now we have our transect xt,yt,zt

    % I want to keep only points with eleveation data
    ffirst = find(~isnan(zt),1,'first');
    flast = find(~isnan(zt),1,'last');

    if ~ isempty(ffirst) && ~isempty(flast)
        xTemp = xt(ffirst:flast);
        yTemp = yt(ffirst:flast);
        zTemp = zt(ffirst:flast);
    else 
        % no data
        continue
    end

    % we need another variable, distance along the transect.
    % we want the transect to begin on land
    % sometimes first or last point is a NaN
    if zTemp(1) > zTemp(end)
        a = 1;
    else
        a = length(zTemp);
    end
    
    % now calcuate distance along transect (probably don't need a loop
    % here, ah well, it works)
    clear dt
    for j = 1 : length(zTemp)
        dt(j) = dist(xTemp(j),yTemp(j),xTemp(a),yTemp(a));
    end

    % I want to sort data by dt
    [dT,I] = sort(dt);
    xT = xTemp(I);
    yT = yTemp(I);
    zT = zTemp(I);
    
    % save the last point on the transect so we an see where it is
    endXYZ(i,1:3) = [xT(end) yT(end) zT(end)];
    
    % now find marsh edge

    % first find (x,y) at MHW
    % do this by finding out which points are on either side of it
    foo = zT - MHW;
    fMHW = find(foo>0,1,'last');
    % MHW is between point #(fMHW-1) and point # fMHW
    % I'll want to get (x,y) at MHW, this will tak a bit of doing
    % first find distance from point fMHW to MHW
    if fMHW > 2
        % tramsect covers MHW
        diffZ = zT(fMHW) - MHW;
        slopeMHWrad = atan((zT(fMHW-1) -zT(fMHW)) / (dT(fMHW-1) - dT(fMHW)));
        dist2MHW = diffZ / sin(slopeMHWrad);
        
        % now back to transect, (use slope of transect in UTM space)
        % used equation for finding (x,y) at a known distance along a known line
        % x = xo + d/sqrt(1+m^2))
        xMHW = xT(fMHW) + dist2MHW/(sqrt(1+perp(i,1).^2));
        % now that know x, use eqn of line
        % y = mx +b
        yMHW = perp(i,1).*xMHW + perp(i,2);
        
        % I want slope in degrees, it makes more sense to me
        slopeMHW = rad2deg(atan((zT(fMHW-1) -zT(fMHW)) / (dT(fMHW-1) - dT(fMHW))));
        
        % check if MHW is too far from profile origin
        if abs(dT(fMHW)) < 100
            % output data
            mhwXYZm(i,1:4) = [xMHW yMHW MHW slopeMHW];
        else
%             keyboard
            disp('possible problem, see line 130')
        end
    else
        % transect does not cover MHW, use slope at top of transect 
        slopeMHW = rad2deg(atan((zT(1) -zT(2)) / (dT(1) - dT(2))));
        fMHW = 1;
    end

    % now see if there is a max in slope seaward of MHW, start by
    % calculating slope
    slope = NaN.*ones(size(zT));
    for j = 1:length(zT)-1
        slope(j) = rad2deg(atan((zT(j+1) -zT(j)) / (dT(j+1) - dT(j))));
    end
  
    % save first point on transect that is below MTL 
    fMTL = find(zT < MTL,1,'first');
    if ~isempty(fMTL)
        mtlXYZm(i,1:4) = [xT(fMTL) yT(fMTL) zT(fMTL) slope(fMTL)];
    end
        
    
    % It seemed easist to me to just cut out the part of the data I needed.
    % I want to stop looking just below MTL
    fMTL = find(zT < (MTL - 0.5),1,'first');
    if isempty(fMTL)
        % no data below MTL
        fMTL = length(zT)-1;
    end
    % OR
%    fMTL = length(zT)-1;
    % start looking at "offset" meters offshore of MHW
    % offset was defined by the user
    ff=find(dT>(dT(fMHW)+offset),1,'first');
    slopeP = slope(ff:fMTL);
    xP = xT(ff:fMTL);
    yP = yT(ff:fMTL);
    zP = zT(ff:fMTL);
    dP = dT(ff:fMTL);
    
    % find minimum slope
    [m,n] = nanmin(slopeP);
    
%     % find any slopes withing 0.5 of this slope
%     temp = slopeP-m;
%     ff = find(temp<2);
%     % use the most seaward one
%     if ~isempty(ff)
% %        n = ff(end);
%         n = floor(mean(ff));
%     end
    
    % check to see if minimum slope is less than slope at MHW 
    % BTW, on the standard forshore, slope will be negative
    if threshold > 0 
        % use wants to use olnly slopes excedding slope at MHW
        if m < (slopeMHW - threshold) 
            % this is a marsh
            xMarsh = xP(n);
            yMarsh = yP(n);
            zMarsh = zP(n);
            slopeMarsh = m;
            % output data
            marshXYZm(i,1:4) = [xMarsh yMarsh zMarsh slopeMarsh];
        else
            % no marsh,  leave as NaN
        end
    elseif ~isempty(m)
        % user wants to use max slope, if one was found
        xMarsh = xP(n);
        yMarsh = yP(n);
        zMarsh = zP(n);
        slopeMarsh = m;
        % output data
        d2r = dist(perp(i,3),perp(i,4),xMarsh,yMarsh);
%        if d2r < 200
        if xMarsh > 1e4
            marshXYZm(i,1:4) = [xMarsh yMarsh zMarsh slopeMarsh];
        else
            keyboard
        end
    else
        % no solution found, leave as NaN
     end

end

%% smooth result

xs = smoothdata(marshXYZm(:,1),'gaussian',5);
ys = smoothdata(marshXYZm(:,2),'gaussian',5);

% some points are bad
dd = dist(marshXYZm(:,1),marshXYZm(:,2),xs,ys);
ff =find(dd > 5000);
xs(ff) = NaN;
ys(ff) = NaN;
% get rid of zeros
ff =find(xs < 1e4);
xs(ff) = NaN;
ys(ff) = NaN;

marshSmoothXY = [xs ys ];

%%

figure
plot(marshXYZm(:,1),marshXYZm(:,2),'r.-')
hold on
plot(xs,ys,'k.-')
legend('raw','smoothed')
title('Marsh shoreline')
axis equal



disp('done')











 