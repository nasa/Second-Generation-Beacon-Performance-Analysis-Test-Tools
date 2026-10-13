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
period=dt;

if length(p0) ~= length(t0) || length(t0) < 2
    error('p0 and t0 must be equal-length vectors containing one complete burst.');
end
if any(diff(t0) <= 0) || max(abs(diff(t0)-1/fs)) > 1e-6/fs
    error('t0 must be monotonically increasing and uniformly sampled.');
end

p0=unwrap(p0);                       % phase must be continuous for differencing
filterSamples=max(3,round(fs*20e-3)); % 20 ms robust phase smoothing
filterSamples=filterSamples+mod(filterSamples+1,2); % use an odd window length
freq=zeros(1,24);
tmeas=zeros(1,24);
tburst=t0-t0(1);

for ix=1:24
    tstart=(ix-1)*period;
    istart=find(tburst >= tstart,1,'first');
    iend=find(tburst >= tstart+dt,1,'first');
    if isempty(istart)
        continue; %error('The input does not contain all 24 frequency measurements.');
    end
    if  isempty(iend)
        iend=length(tburst);
    end
    elapsed=t0(iend)-t0(istart);
    phaseWindow=p0(istart:iend);
    %phaseWindow=movmedian(phaseWindow,filterSamples,'Endpoints','shrink');
    freq(ix)=(phaseWindow(end)-phaseWindow(1))/elapsed/2/pi;
    tmeas(ix)=(tburst(istart)+tburst(iend))/2;
end

d2=zeros(1,20);
for ix=1:20
    d2(ix)=max(freq(ix:ix+4))-min(freq(ix:ix+4));
end

sfigure(fig);clf; hold on;
offset=mean(freq);
plot(tmeas,freq-offset,'-x');
tx=((1:20)+1)*period+tmeas(1);
[y,idx]=max(d2);
plot(tx,d2,'-x');
plot(tx(idx),y,'go');
grid on;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('Short Frequency Variation=%0.3fppb', y/406.050e6*1e9));
legend('Frequency', 'Max-Min', 'Worst','location', 'best');

dev=y/406.050e6*1e9; %ppb

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
rng=fix(0.01*fs:fix(period*fs):0.99*fs);
plot(t0(rng(2:end)),diff(res(rng))./diff(t0(rng))/2/pi,'-x');title('Freq Variation vs Time');xlabel('Time (s)'); ylabel('Frequency (Hz)'); grid on;

