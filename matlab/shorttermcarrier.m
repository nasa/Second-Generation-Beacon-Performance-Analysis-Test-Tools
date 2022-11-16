% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: shorttermcarrier 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 12 $ $Date: 2020-01-21 08:10:30 -0500 (Tue, 21 Jan 2020) $
% ///				$Id: shorttermcarrier.m 12 2020-01-21 13:10:30Z reesebo $
% ///            
function [dev,offset,d2]=shorttermcarrier(p0,t0,dt,fig)

fs=1/(t0(2)-t0(1));
%dt=166.66667e-3;
fx=fix(fs*dt);


p0=medfilt1(p0,fx,10);

win=[1 fix(fx)]+fix(fx/2);

try
for(ix=1:length(p0))
    dev(ix) = diff(p0(win))/dt/2/pi;
    win=win+fx;
end
catch
    %dev(ix) = diff(p0([win(1) length(p0)-fix(fx/10)]))/dt;
end


try
win=1:5;
for(ix=1:length(dev))
    m=max(dev(win));
    n=min(dev(win));
    d2(ix)=m-n;
    win=win+1;
end
catch
end

sfigure(fig);clf; hold on;
tx=(1:length(dev))*dt;
offset=mean(dev);
plot(tx,dev-offset,'-x');
tx=(1:length(d2))*dt+5*dt/2;
[y,idx]=max(d2);
plot(tx(idx),y,'go');
plot(tx,d2,'-x');
grid on;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('Short Frequency Variation=%0.3fppb', y/406*1e3));
legend('Frequency', 'Max-Min', 'Worst','location', 'best');

dev=y/406*1e3; %ppb

sfigure(fig+900);clf;
subplot(311);hold on;
rng=fix(0.01*fs:0.99*fs);
plot(t0(rng),p0(rng));title('Phase vs Time');xlabel('Time (s)'); ylabel('Phase (rads)');grid on;
p=polyfit(t0(rng),p0(rng),1); % note linear trend in phase is fixed frequency.
plot(t0(rng),polyval(p,t0(rng)));
subplot(312);
res=p0-polyval(p,t0);
plot(t0(rng),res(rng));title('Phase Variation vs Time (residual)');xlabel('Time (s)'); ylabel('Phase (Rads)');grid on;
subplot(313);
rng=fix(0.01*fs:fix(41e-3*fs):0.99*fs);
plot(t0(rng(2:end)),diff(res(rng))/dt/2/pi,'-x');title('Freq Variation vs Time');xlabel('Time (s)'); ylabel('Frequency (Hz)'); grid on;

