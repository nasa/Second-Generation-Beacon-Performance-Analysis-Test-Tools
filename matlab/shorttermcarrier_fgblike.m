% ///  Copyright(c) 2017 by Reese Bovard, all rights reserved
% ///  Proprietary Information of Concentric Real Time, LLC
% ///  Subject to Non-Disclosure.  Do Not Copy or Distribute
% ///  
% ///     Project: T21 Software Tools for Cospas-Sarsat DSSS Beacons
% ///     
% ///     Module: shorttermcarrier 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision$ $Date$
% ///				$Id$
% ///            
function [dev,s2,d2]=shorttermcarrier_fgblike(p0,t0,fig)

fs=1/(t0(2)-t0(1));
dt=100e-3;
fx=fix(fs*dt);


p0=medfilt1(p0,fx);

win=[1 fix(fx)]+fix(fx/2);

try
for(ix=1:length(p0))
    dev(ix) = diff(p0(win))/dt;
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
s2=dev(2);
dev=abs(dev(4)-dev(2))/406*1e3; %ppb

title(sprintf('Short Frequency Variation=%0.3fppb', dev));
legend('Frequency', 'Max-Min', 'Worst','location', 'best');

%dev=dev(
sfigure(fig+900);clf;
subplot(311);hold on;
rng=fix(0.01*fs:0.99*fs);
plot(t0(rng),2*pi*p0(rng));title('Phase vs Time');xlabel('Time (s)'); ylabel('Phase (rads)');grid on;
p=polyfit(t0(rng),p0(rng),1); % note linear trend in phase is fixed frequency.
plot(t0(rng),2*pi*polyval(p,t0(rng)));
subplot(312);
res=p0-polyval(p,t0);
plot(t0(rng),2*pi*res(rng));title('Phase Variation vs Time (residual)');xlabel('Time (s)'); ylabel('Phase (Rads)');grid on;
subplot(313);
rng=fix(0.01*fs:fix(41e-3*fs):0.99*fs);
plot(t0(rng(2:end)),diff(res(rng))/dt,'-x');title('Freq Variation vs Time');xlabel('Time (s)'); ylabel('Frequency (Hz)'); grid on;

