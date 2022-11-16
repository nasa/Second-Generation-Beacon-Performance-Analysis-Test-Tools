% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: burstdetection 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: burstdetection.m 15 2022-09-29 15:45:13Z reesebo $
% ///            


function t0=burstdetection(t,p,fig)


T=(t(2)-t(1));
fs=fix(1/T);
twin=1e-3;

winsz=round(twin*fs);  % 1ms window
win=1:winsz;
ts=twin/2;
for(ix=1:length(p))
    try
    ox(ix) = max(p(win));
    om(ix) = min(p(win));
    oa(ix) = mean(p(win));
    tx(ix) = ts+(ix-1)*T;
    win=win+1;
    catch
        break;
    end
end

sfigure(fig);clf; hold on;

plot(tx,ox);
plot(tx,om);
plot(tx,oa);

km=kmeans(om(:),2);

idx1=find(km==1);  % find cluster of greater value
idx2=find(km==2);

if(mean(om(idx1)) > mean(om(idx2)))
    plot(tx(idx1),om(idx1),'go');
    t0=tx(idx1(1)); %-twin/2;;
    mx=max(ox(idx1));
    mn=min(om);
    plot([t0 t0], [mn mx],'r','LineWidth',4);
else
    plot(tx(idx2),om(idx2),'go');
    t0=tx(idx2(1)); %-twin/2;
    mx=max(ox(idx2));
    mn=min(om);
    plot([t0 t0], [mn mx],'r','LineWidth',4);
end
title('Power over time/ burst detection');
xlabel('time (s)');
ylabel('power (dB)');

