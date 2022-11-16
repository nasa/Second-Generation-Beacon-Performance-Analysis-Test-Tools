% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: doevm 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: doevm.m 15 2022-09-29 15:45:13Z reesebo $
% ///            
function [mxx,ta,ua]=doevm(fchip,rx,fig)
rx=rx(:);

%winsz=round(150e-3*fchip);
winsz=round(5e-3*fchip);
stepsz=round(1e-3*fchip);
t0 = stepsz/fchip;
strt = fix(1e-3*fchip);
endx = fix(fchip-winsz-strt);
win=1:winsz;
%rx=[rx(1:end-strt) ;rx(win+endx-1)];
win = win+strt;

con =[   0.7071 + 0.7071i; ...
  -0.7071 + 0.7071i; ...
  -0.7071 - 0.7071i; ...
   0.7071 - 0.7071i];

evm = comm.EVM( ...
    'ReferenceSignalSource','Estimated from reference constellation', ...
    'ReferenceConstellation',con, ...
    'MeasurementIntervalSource', 'Custom with periodic reset','MeasurementInterval',winsz);
mx=[];t=[];
for(ix=1:2*fchip/winsz*1.5*10)
try
    mx(ix)=step(evm,rx(win));
    t(ix)=ix*t0;
    win=win+stepsz;
catch
    break;
end

end
% find signal acquisition in 150ms
try
ax = find(mx(1:150)>15);
ua = mean(mx(1:150));
catch
    ax=[];
    ua=0;
end
if(~isempty(ax))
    ta=t(ax(end));
else
    ta=0;
end
idx = find(t>0.01 & t< 0.99);
mx=mx(idx);
t=t(idx);
try
    win=1:150;
for(ix=1:length(mx))
   mxx(ix) = mean(mx(win)); 
   tx(ix) = 75e-3+ix*10e-3;
   win=win+10;
end
catch
    
end
sfigure(fig);clf;hold on;
plot(t,mx);
plot(tx,mxx,'ro-');
grid on;
legend('EVM%', 'uEVM%');
title('EVM over time');
xlabel('time (s)');
ylabel('EVM %');