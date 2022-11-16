% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: getwfPN23_upsamplewithdata 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: getwfPN23_upsamplewithdata.m 15 2022-09-29 15:45:13Z reesebo $
% ///            
% makewfPN23 Oversampled/Filtered.
function [xm, xs]=getwfPN23_upsamplewithdata(nSamples, vsg_filename, pulsetype, mode, data)
% mode='selftest' or 'normal';

%% build waveform NRZ with lookup table using sinc pulse

fchip=38400;
PL = 6400;
N  = 250;

[ si, sq ] = getPN23( fchip, mode );

% apply modulation to sequence
siq = zeros(1,2*fchip+2);
siq(1:2:end-3) = si;
siq(2:2:end-2) = sq;
siq(end-1) = si(end);
siq(end)   = sq(end);

miq = siq;
r = [ 2*PL+1 : 2 : 2*PL+512 ];
for n=0:N/2-1
    miq(r  ) = xor( miq(r  ), data(2*n+1) ); 
    miq(r+1) = xor( miq(r+1), data(2*n+2) ); 
    r = r + 512;
end;

% modulated and spread
I=2*[miq(1:2:end)]-1;  % data NRZ
Q=2*[miq(2:2:end)]-1;

xm=do_oqpsk_upsample(I,Q,nSamples,fchip,pulsetype);

% spread sequence
I=2*[siq(1:2:end)]-1;  % data NRZ
Q=2*[siq(2:2:end)]-1;

xs=do_oqpsk_upsample(I,Q,nSamples,fchip, pulsetype);

if(~isempty(vsg_filename))
    make_vsg(xm,vsg_filename);
end

end

function oo=do_oqpsk_upsample(I,Q,nSamples,fchip,pulsetype)
doplot=false;
ncycles=.5;

rc=sqrt(rcosdesign(0.25,1,nSamples)); % narrower, higher amplitude, sharper IQ corners

rc=rc/max(rc);
rc=rc(1:end-1);

switch(pulsetype)
    case 'sin'
        LUT_pulse=(sin( (0:nSamples-1) * 2*pi*0.5/nSamples));

    case 'rrc'
        LUT_pulse=rc;
        
    case 'tri'
        LUT_pulse= [[1:nSamples/2] [nSamples/2+1:-1:2]]/nSamples;
        
    case 'none'
        LUT_pulse = 0.7071*ones(1,nSamples);

end
win=1:nSamples;
for(ix=1:length(I))
    oi(win) = I(ix) * LUT_pulse;
    oq(win) = Q(ix) * LUT_pulse;
    win=win+nSamples;
end
oo=[oi(:); zeros(1,nSamples/2)']+1j*[zeros(1,nSamples/2)' ; oq(:)];
oo=oo.';

if(doplot)
    sfigure(3);clf;plot(LUT_pulse);
    title('Pulse Shape');
    sfigure(5);
    plot(oo(nSamples:end-nSamples),'.-')
    title('IQ');
    sfigure(7);
    quickfft(oo,fchip*nSamples,'b');
end
end