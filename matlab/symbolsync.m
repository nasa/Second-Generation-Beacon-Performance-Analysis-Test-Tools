% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: symbolsync 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: symbolsync.m 15 2022-09-29 15:45:13Z reesebo $
% ///            
function [rxSync,err]=symbolsync(y,fs,sps,tp,fig)
%y=y/max(abs(y));

% rot=unwrap(angle(mean(y(1:4*50))));
% y=y.*exp(1j*rot);

switch(tp)
    case 1
        type={'TimingErrorDetector', 'Zero-Crossing (decision-directed)'};
    case 2
        type={'TimingErrorDetector', 'Gardner (non-data-aided)'};
    case 3
        type={'TimingErrorDetector', 'Early-Late (non-data-aided)'};
    case 4
        type={'TimingErrorDetector', 'Mueller-Muller (decision-directed)'};
end

symbolSync = comm.SymbolSynchronizer('SamplesPerSymbol',sps, 'Modulation', 'OQPSK', type{1}, type{2});
% symbolSync.DampingFactor=1.2;
% symbolSync.NormalizedLoopBandwidth=0.02;
y=y/max(abs(y));
% rr=isnan(y);
% y(rr)=0;
[rxSync,err] = step(symbolSync,y(:));

sfigure(1100);
plot(err/4);
xlabel('Samples');
ylabel('Symbol Error');
% p=polyfit(tx,fx,1);

%% adjust amplitude
winsz=sps;
win=1:winsz;
lx=length(rxSync);
rxSync=[rxSync ;ones(winsz,1)*rxSync(end)];

for(ix=1:lx)
    rxSync1(ix) = rxSync(ix)/mean(abs(rxSync(win)));
    win=win+1;
end
rxSync=rxSync1;
if(fig)
    sfigure(fig);clf;hold on;
    
    %plot(rxSync(fix(length(rxSync)/2):end)); hold on;
    plot(rxSync(fix(length(rxSync)/2):end),'.');
    Constellation=.7071*[1+1i,1-1i,-1+1i,-1-1i];

    plot(Constellation,'gd');
    title('symbol map center to end');
    axis([-1 1 -1 1]);
    grid;

    sfigure(fig+1); clf; hold on;

    %plot(rxSync); hold on;
    plot(rxSync,'.'); 
    plot(Constellation,'gd');
    title('symbol map all');
    axis([-1 1 -1 1]);
    grid;
end