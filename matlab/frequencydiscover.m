% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: frequencydiscover 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: frequencydiscover.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

function [t,f,p]=frequencydiscover(xo,fs,offset,pcal,fig)



%stepsz=round(65e-6*fs);
stepsz=round(1e-3*fs);
n = 50*stepsz;
N = 2.^(fix(log10(n)/log10(2)));
win=1:N;
ferr=[];
lwin=length(win);
for(ix=1:length(xo)/stepsz-1)
    win=win+stepsz;
    try
        [centerfreq, ~, CN0, CPow, ~]=specpowr_sr2(xo(win),-10000,10000,fs,pcal,0);
    catch e
        %disp(e);
        break;
    end
    t=win(lwin/2)/fs;
    ferr(ix,:) = [t centerfreq/4 real(CN0)/4 CPow]; % note cno is computed with 20log10 in specpeak_sr2, so we have an extra power of 4 to remove
    if(mod(ix,1000)==0)
        fprintf('.');
    end
end

t=ferr(:,1);
f=ferr(:,2);
%p=ferr(:,3);
p=ferr(:,4);

if(0)
for(ix=1:10)
    u=mean(f);
    s=std(f);
    idx=find(abs(f)<abs(u+6*s));
    t=t(idx);
    p=p(idx);
    f=f(idx);
end
end

if(fig)
    sfigure(fig);clf;
    subplot(211);
    plot(t,f+offset);title('freq over time');
    xlabel('time (s)');
    ylabel('freq (Hz)');
    subplot(212);
    plot(t,p);title('power over time');
    xlabel('time (s)');
    ylabel('power dBm');
    drawnow;
end