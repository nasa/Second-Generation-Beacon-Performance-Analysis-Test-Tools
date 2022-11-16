% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: getwfPN23_upsample 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: getwfPN23_upsample.m 15 2022-09-29 15:45:13Z reesebo $
% ///            
% makewfPN23 Oversampled/Filtered.
function [xm, xs]=getwfPN23_upsample(nSamples, vsg_filename, pulsetype, mode)
% mode='selftest' or 'normal';

%% build waveform NRZ with lookup table using sinc pulse

fchip=38400;
PL = 6400;
N  = 250;

[ si, sq ] = getPN23( fchip, mode );

% % set data
canned={ '01234567', '89abcdef', 'deadbeef', '01234567', '89abcdef', 'deadbeef', 'FFFc2c15', '9136dcc0'};

data=[];
for(ix=1:length(canned))
    data=[data dec2bin(hex2dec(canned{ix}),32)];
end
data=[str2bin(data(1:250))];

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

xm=do_oqpsk_upsample(I,Q,nSamples,fchip,pulsetype,0);

% spread sequence
I=2*[siq(1:2:end)]-1;  % data NRZ
Q=2*[siq(2:2:end)]-1;

xs=do_oqpsk_upsample(I,Q,nSamples,fchip, pulsetype,0);

if(~isempty(vsg_filename))
    make_vsg(xm,vsg_filename);
end

end

function oo=do_oqpsk_upsample(I,Q,nSamples,fchip,pulsetype, doplot)




switch(pulsetype)
    case 'msk'
        LUT_pulseup=0.7071*(sin( (0:nSamples-1) * 2*pi*0.5/nSamples));
        LUT_pulsedn=-0.7071*(sin( (0:nSamples-1) * 2*pi*0.5/nSamples));
    case 'cos'
        sind=(cos( (0:nSamples-1) * 2*pi/nSamples/2))*0.7071;
        LUT_pulsedn=sind;
        LUT_pulseup=-1*LUT_pulsedn;
        sc=0.7071;
    case {'rrc','rrcp2'}    
        mf=mfp2(fchip*nSamples,nSamples-2);
        rc=mf.Numerator;
        sig=[ones(1,nSamples) -1*ones(1,nSamples)];
        op=filter(rc,1,sig);
        % figure(4);clf;hold on;
        % sc=1/max(op);
        % op=op*sc;
        % plot(t,op,'x-');
        % plot(t,sig+8/fs);
        idx=nSamples:2*nSamples-1;
        % plot(t(idx),op(idx),'o');

        rcdn=op(idx)*1/max(op(idx))*0.7071;
        rcup=fliplr(rcdn);

        LUT_pulseup=rcup;
        LUT_pulsedn=rcdn;
        sc=0.7071;
    case 'rrcp8'    
        mf=mfp8(fchip*nSamples,nSamples-2);
        rc=mf.Numerator;
        sig=[ones(1,nSamples) -1*ones(1,nSamples)];
        op=filter(rc,1,sig);
        % figure(4);clf;hold on;
        % sc=1/max(op);
        % op=op*sc;
        % plot(t,op,'x-');
        % plot(t,sig+8/fs);
        idx=nSamples:2*nSamples-1;
        % plot(t(idx),op(idx),'o');

        rcdn=op(idx)*1/max(op(idx))*0.7071;
        rcup=fliplr(rcdn);

        LUT_pulseup=rcup;
        LUT_pulsedn=rcdn;
        sc=0.7071;        
    case 'tri'
        LUT_pulsedn= linspace(1,-1,nSamples)*0.7071;
        LUT_pulseup= -1*(LUT_pulsedn);
        sc=0.7071;
    case 'rec'
        LUT_pulseup = 0.7071*ones(1,nSamples);
        LUT_pulsedn = -0.7071*ones(1,nSamples);
        sc=0.7071;
    case {'bessel'}
        
        fs=fchip*nSamples;
        w0=50000;
        w=w0/fs;

        [b,a]=besself(nSamples-2,w);

        sig=[ones(1,nSamples) -1*ones(1,nSamples)];
        op=filter(b,a,sig);
        % figure(4);clf;hold on;
        % sc=1/max(op);
        % op=op*sc;
        % plot(t,op,'x-');
        % plot(t,sig+8/fs);
        idx=nSamples+1:2*nSamples;
        % plot(t(idx),op(idx),'o');

        rcdn=op(idx)*1/max(op(idx))*0.7071;
        rcup=fliplr(rcdn);

        LUT_pulseup=rcup;
        LUT_pulsedn=rcdn;
        sc=0.7071;
    case 'none'
        LUT_pulseup = ones(1,nSamples);
        LUT_pulsedn = -1*ones(1,nSamples);
        sc=1;


end
win=1:nSamples;
Is=0;
Qs=0;

if(strcmp(pulsetype, 'msk') == 0)
    for(ix=1:length(I))
        dI=I(ix)-Is;
        dQ=Q(ix)-Qs;
        switch(sign(dI))
            case -1
                oi(win) =LUT_pulsedn;
            case 1
                oi(win) =LUT_pulseup;
            case 0
                oi(win) = sc*I(ix);
        end
        switch(sign(dQ))
            case -1
                oq(win) =LUT_pulsedn;
            case 1
                oq(win) =LUT_pulseup;
            case 0
                oq(win) =sc*Q(ix);
        end
        Is=I(ix);
        Qs=Q(ix);
        win=win+nSamples;
    end
else
    for(ix=1:length(I))
        dI=I(ix);
        dQ=Q(ix);
        switch(sign(dI))
            case -1
                oi(win) =LUT_pulsedn;
            case 1
                oi(win) =LUT_pulseup;
        end
        switch(sign(dQ))
            case -1
                oq(win) =LUT_pulsedn;
            case 1
                oq(win) =LUT_pulseup;

        end

        win=win+nSamples;
    end
end
oo=[oi(:); zeros(1,nSamples/2)']+1j*[zeros(1,nSamples/2)' ; oq(:)];
oo=oo.';
ob=obw(oo,nSamples*fchip);
fprintf('OBW = %0.4f\n', ob);
ob=ob/1e3;
if(doplot)
    sfigure(3);clf;plot(LUT_pulseup);
    title('Pulse Shape');
    sfigure(5);
    plot(oo(nSamples:end-nSamples),'.-')
    title('IQ');
    sfigure(7);
    quickfft_z_100(oo,fchip*nSamples);
    patch([-ob/2 ob/2 ob/2 -ob/2], [0, 0, -100 -100],'r');alpha(0.3);
    title(sprintf('%s: occupied BW=%dkHz', pulsetype, fix(ob)));
end
end