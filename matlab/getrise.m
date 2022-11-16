% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: getrise 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 12 $ $Date: 2020-01-21 08:10:30 -0500 (Tue, 21 Jan 2020) $
% ///				$Id: getrise.m 12 2020-01-21 13:10:30Z reesebo $
% ///            
function [rt, t90, t10]=getrise(xb,fs,fig)

l2 = 0.5*fs;
p=abs(xb(:));
t=(0:length(xb)-1)/fs;
%p=medfilt1(p,3);

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
plot(t(idx1),p(idx1),'x');
plot(t(idx2),p(idx2),'.');
m1=mean(p(idx1));
m2=mean(p(idx2));

if(m1>m2)
    mx=m1;
    mn=m2;
    idxmx = idx1;
    idxmn = idx2(idx2<l2);
else
    mx=m2;
    mn=m1;
    idxmx = idx2;
    idxmn = idx1(idx1<l2);
end

thresh=(mx-mn)/2;
 plot([t(idxmx(1)) t(idxmx(end))], [mx mx]);
 plot([t(idxmn(1)) t(idxmn(end))], [mn mn]);
 plot([t(1) t(end)], [thresh thresh]);
 legend('Signal','Absent','Present','uMax', 'uMin', 'Threshold')
 xlabel('Time (s)')
 ylabel('Signal Voltage');

step= (thresh-mn)/1000;
 while(idxmx(1) < idxmn(end)) % overlap
    thresh=thresh-step;
    idxmx = find(p(t<.5)>thresh);
    idxmn = find(p(t<.5)<=thresh);
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

rise = [idxmn(end-mns:end); idxmx(1:samps)];

sfigure(fig);clf;hold on;
plot(t(rise), p(rise),'x-');
plot([t(rise(1)) t(rise(end))], [pn pn]);
plot([t(rise(1)) t(rise(end))], [thmx thmx]);
plot([t(rise(1)) t(rise(end))], [thmn thmn]);

brk = find(p(rise)>thmn);
plot(t(rise(brk(1))), p(rise(brk(1))),'o');
plot(t(rise(brk(1)-1)), p(rise(brk(1)-1)),'o');
px=polyfit([p(rise(brk(1))) p(rise(brk(1)-1))],[t(rise(brk(1))) t(rise(brk(1)-1))],1);
t0=polyval(px,thmn);
t10=t0;

brk = find(p(rise)>thmx);
plot(t(rise(brk(1))), p(rise(brk(1))),'o');
plot(t(rise(brk(1)-1)), p(rise(brk(1)-1)),'o');
px=polyfit([p(rise(brk(1))) p(rise(brk(1)-1))],[t(rise(brk(1))) t(rise(brk(1)-1))],1);
t1=polyval(px,thmx);
rt=t1-t0;

t90=t1;
title(sprintf('Rise Time = %0.2g', rt));
ylabel('Signal Amplitude');
xlabel('Time(s)');
legend('Signal', 'Pn', '0.9*Pn', '0.1*Pn')

complete=true;
catch
    iter=iter-1;
    continue;
end
end