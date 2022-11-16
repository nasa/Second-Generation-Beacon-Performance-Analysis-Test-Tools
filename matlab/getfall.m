% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: getfall
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: getfall.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
function [rt, t90]=getfall(xb,fs,fig)
l2 = .2*fs; %fix(length(xb)/2);

p=abs(xb(end-l2:end));p=p(:);
t0=(0:length(xb)-1)/fs;
t=t0(end-l2:end);
complete=false;
iter=10;
rt=0;
t90=0;

while(~complete && iter>0)

km=kmeans(p,2,'Replicates',10);

try
idx1=find(km==1);  % find cluster of greater value
idx2=find(km==2);

sfigure(fig); hold on;
plot(t,p);
plot(t(idx1),p(idx1),'.');
plot(t(idx2),p(idx2),'.');
m1=mean(p(idx1));
m2=mean(p(idx2));

if(m1>m2)
    mx=m1;
    mn=m2;
    idxmx = idx1;%(idx1>l2);
    idxmn = idx2;%(idx2>l2);
else
    mx=m2;
    mn=m1;
    idxmx = idx2;%(idx2>l2);
    idxmn = idx1;%(idx1>l2);
end
thresh=(mx-mn)/2;
 plot([t(idxmx(1)) t(idxmx(end))], [mx mx]);
 plot([t(idxmn(1)) t(idxmn(end))], [mn mn]);
 plot([t(1) t(end)], [thresh thresh]);
step= (thresh-mn)/1000;
 while(idxmx(end) > idxmn(1)) % overlap
    thresh=thresh-step;
    idxmx = find(p(t>.5)>thresh);%+l2;
    idxmn = find(p(t>.5)<=thresh);%+l2;
 end
 
pn=mx;

thmn=0.1*pn; % min thresh
thmx=0.9*pn; % max thresh

samps = fix(5e-4*fs);

if(samps > length(idxmn))
    mns = length(idxmn)-1;
else
    mns = samps;
end

fall = [idxmx(end-samps: end); idxmn(1:mns)];

sfigure(fig);clf;hold on;
plot(t(fall), p(fall),'x-');
plot([t(fall(1)) t(fall(end))], [pn pn]);
plot([t(fall(1)) t(fall(end))], [thmx thmx]);
plot([t(fall(1)) t(fall(end))], [thmn thmn]);

brk = find(p(fall)<thmn);
plot(t(fall(brk(1))), p(fall(brk(1))),'o');
plot(t(fall(brk(1)-1)), p(fall(brk(1)-1)),'o');
px=polyfit([p(fall(brk(1))) p(fall(brk(1)-1))],[t(fall(brk(1))) t(fall(brk(1)-1))],1);
t0=polyval(px,thmn);

brk = find(p(fall)>thmx);
plot(t(fall(brk(end))), p(fall(brk(end))),'o');
plot(t(fall(brk(end)+1)), p(fall(brk(end)+1)),'o');
px=polyfit([p(fall(brk(end))) p(fall(brk(end)+1))],[t(fall(brk(end))) t(fall(brk(end)+1))],1);
t1=polyval(px,thmx);
rt=t0-t1;

t90=t1;
title(sprintf('Fall Time = %0.2g', rt));
ylabel('Signal Amplitude');
xlabel('Time(s)');
legend('Signal', 'Pn', '0.9*Pn', '0.1*Pn')

complete=true;
catch e
    getReport(e)
    iter=iter-10;
    continue;
end
end