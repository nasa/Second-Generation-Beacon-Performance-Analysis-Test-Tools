% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: getpayload 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 16 $ $Date: 2022-10-14 10:17:46 -0400 (Fri, 14 Oct 2022) $
% ///				$Id: getpayload.m 16 2022-10-14 14:17:46Z reesebo $
% ///            
function specs = getpayload(x,fsi,risetime,pcal)

specs.error=false;
specs.errormessage='';

modes={'normal', 'selftest'};
mftypes={'none', 'rrc', 'sin', 'sinx'};

fchip=38400;

fig=1;
mx0=inf;
ta0=inf;
ua0=inf;
eps0=0.5;
mxx=[];
rxx=[];
evmb=zeros(1,3);
snrb=zeros(1,3);
evmbr=zeros(1,3);
snrbr=zeros(1,3);

pkamp=zeros(1,3);
for(iiz=1:3)
    bn0=0.1;
    mf_type=mftypes{iiz};        

    switch(mf_type)
        case 'none'
            fs=153600;
            disp('No matched filter on');
        case 'rrc'
            fs=153600;
            disp('rrc matched filter on');
            x=matchfilter(x,fsi,fchip);
        case 'sin'
            fs=153600;
            x=halfsinepulse(x,fsi,fchip);
            bn0=1;
            disp('msk processing');      

    end

fprintf('profile power and frequency...');
[t,f,p]=frequencydiscover(x,fsi,0,pcal,fig); 
saveasrsb(fig,'powerfreq.png');fig=fig+1;
fprintf('complete\n');
figpwr=fig-1;
fprintf('detect burst energy...');
t0=burstdetection(t,p,fig); 
saveasrsb(fig,'burstdetect.png');fig=fig+1;
fprintf('complete\n');
t0=t0-risetime;

t0s=fix(t0*fsi);
if(t0s<0)
    t0s=1;
    t0=0;
end

t0send=t0s+fix(fsi*1.2)-1;
if(t0send> length(x))
    t0send=length(x);   
end
xb=x(t0s:t0send); 

if(t0send-t0s<fsi)
    specs.error=true;
    specs.errormessage='burst length too short';
    return;
end

fprintf('resample...');

sps=fs/fchip;
[a,b]=rat(fsi/fs);
xob=resample(xb,b,a);
fprintf('complete\n');

fprintf('get rise/fall time...');
try
[rt,ts90,ts10]=getrise(xob,fs,fig); 
saveasrsb(fig,'risetime.png');fig=fig+1;
[ft,te90]=getfall(xob,fs,fig); 
saveasrsb(fig,'falltime.png');fig=fig+1;
fprintf('complete\n');
catch
specs.error=true;
specs.errormessage='testing rise';
return;
end


dur=te90-ts90;
%dur=1;
if(dur < 0.99)
    specs.error=true;
    specs.errormessage='burst length too short';
    return;
end

% capture power info
sfigure(figpwr);
subplot(211);
axis([ts90 te90+0.1 -1000 1000]);
subplot(212);
axis([ts90 te90+0.1 26 27]);
idx=find(t>ts90+.1 & t<te90);
u=mean(p(idx));
s=std(p(idx));
mx=max(p(idx));
mn=min(p(idx));
specs.power.u=u;
specs.power.mx=mx;
specs.power.mn=mn;
specs.power.s=s;
drawnow;



st=max([1,fix((ts10-risetime)*fs)]);
if(st>1)
    t0 = t0+st/fs;
end

xob=xob(st:end);
fprintf('profile power and frequency over burst...');
[xoc,cfx,cr]=center_scale(xob,fchip,fs,fig);
saveasrsb(fig,'centerscale.png');fig=fig+1;
fprintf('complete\n');
xoc=xoc/max(abs(xoc));


res=[];
for(ix=1:length(modes))
    mode=modes{ix};

    [res,t1,pk]=getdetect(xoc,fs,fchip,mode);
    if(~isempty(res))
        xoc=res;
        break;   
    end
end
if(isempty(res))
    specs.error=true;
    specs.errormessage = 'no signal detected';
    return;
end
pkamp(iiz) = pk;
t0=t0+t1;

fs=2*fs;
xbb=x(round(t0*fsi):end).';
tt=0:length(xbb)-1;
xoc = xbb.*exp(-cfx/fsi*2*pi*1j*tt);
sps=fs/fchip;
[a,b]=rat(fsi/fs);
xoc=resample(xoc,b,a);
xoc=xoc(1:fs+2*sps);

y=xoc;

fprintf('carrier phase correction...');
bn=0.0001;
fprintf('carrier phase correction...');
[y,p00,t00]=carrierphasecorrection(y,fs,bn,fchip,fig);
p0x=p00;
t0x=t00;

saveasrsb(fig,'cpc.png');fig=fig+1;
fprintf('complete\n');
% 




[~, xs]=getwfPN23_upsample(sps,'','none', mode); % get reference waveform

fprintf('demod...');
try
    [bits,bsyms,rot,pyy]=dodemodds(y(1:end),xs,sps,fig);
    [cno,evmb(iiz),snrb(iiz)] = docno(bsyms);
    fprintf('bits: evm=%0.2f%%, cno=%0.2f, snr=%0.2f\n', evmb(iiz)*100, cno, snrb(iiz));

    bin2hex([bits(51:end) ],1)
    saveasrsb(fig,'demod.png');fig=fig+1;
catch e
    bits=-1;
    specs.error=true;
    specs.errormessage='signal could not be successfully demodulated';
    return;
end
fprintf('complete\n');

%MSK: Try removing the 1/2 Tb offset before processing.
        
    for(il=1:3)
        y=xoc*exp(1j*rot);
        bn=(0.01)^il*bn0;

        fig=100;
        fprintf('carrier phase correction...');
        [y,p0x,t0x]=carrierphasecorrection(y,fs,bn,fchip,fig);fig=fig+1;
        fprintf('complete\n');

        fprintf('symbolsync...');
        [rx,err]=symbolsync(y,fs,sps,4,fig);fig=fig+1;
        fprintf('complete\n');

        fprintf('correct tilt...');
        rx=correcttilt(rx,fig);fig=fig+1;
        fprintf('complete\n');

        fprintf('calculate EVM...');
        [mx,ta,ua]=doevm(fchip,rx(~isnan(rx)) ,fig); fig=fig+1;
        fprintf('complete\n');
        
        if(ta<=ta0 && ua<ua0 && max(mx)+eps0 < mx0)
            mx0=max(mx);
            ta0=0.075;
            ua0=ua;
            mxx=mx;
            rxx=rx;
            xocc=y;
            bnxx=bn;
            besterr=err;
            using_mf = mf_type;
            p00=p0x;
            t00=t0x;
            saveasrsb(100,'cpc.png');
            saveasrsb(101,'symbsync.png');
            saveasrsb(102,'tilt.png');
            saveasrsb(103,'EVM.png');
            saveasrsb(1100,'Tdchiperr.png');
        end
    end
end

[cno,evmo,snr] = docno(rxx);
fprintf('chips: evm=%0.2f%%, cno=%0.2f, snr=%0.2f\n', evmo*100, cno, snr);

fprintf('chosen bn=%04f\n', bnxx);
if(mod(length(rxx),2) == 0)
    rxx=[rxx 0];
    disp('Added a zero to rx, not divisible by 2!');
end

fprintf('demod after symbol sync\n');
    [ si, sq ] = getPN23( fchip, mode );
    xs = (2*si-1) + 1j* (2*sq-1);

    evmb_ss=inf;
    for(ix=1:3)
        [~,bsyms,~,~]=dodemodds(rxx,xs,1,fig,ix-1);fig=fig+1;
        [mx]=evmsimple(bsyms(~isnan(bsyms)));
        if(mx<evmb_ss)
            evmb_ss=mx;
            rxsav=rx;
        end
    end
fprintf('PN Code Demod...');
for(ixx=1:2)
    pn=qpskdemod(rxx);
    pnb=str2bin(getpnbin(pn));
    [pni,pnq,invi,invq,pki,pkq,posi,posq]=assess_pn(bits,pnb,fchip,mode);

    if(pki>38300 && pkq>38300)
        fprintf('I error positions: [%d]\n', posi);
        fprintf('Q error positions: [%d]\n', posq);
        break;
    else
        rx=rx(2:end-1);
    end

end
fprintf('complete\n');fig=fig+4;

fprintf('chip rate scan...');
tw=0.01; %/fchip*5000; %/2; %10e-3; %
crwin=fix(tw*fs);
step=fix(1/fchip*100*fs);
[crvec, dv]=evalchiperror(besterr,fs,fchip,tw,fig); 
saveasrsb(fig,'ChipError.png');fig=fig+2;%
tw=1/fchip*5000; %/2; %10e-3; %
crwin=fix(tw*fs);
step=fix(1/fchip*100*fs);
[~, dv2]=chipscan(xoc,fs,fchip,crwin,step,fig);saveasrsb(gcf,'chipscan.png');fig=fig+1;

p2p = pk2pk(xocc,fs,fchip,fig);fig=fig+1;
fprintf('complete\n');

[dev_100,s2,maxMinusMin]=shorttermcarrier_fgblike(p00,t00,fig);fig=fig+1;
[dev_167,o2,maxMinusMin]=shorttermcarrier(p00,t00,167.7e-3,fig);
[dev_42,o2,maxMinusMin]=shorttermcarrier(p00,t00,41.66667e-3,fig);
saveasrsb(fig,'shorttermcarrier.png');fig=fig+1;
fprintf('complete\n');


specs.sq.offset=cfx+o2;
specs.sq.dev100=dev_100;
specs.sq.dev167=dev_167;
specs.sq.dev=dev_42;
specs.sq.s2=s2;
specs.sq.maxMinusMin = maxMinusMin;

%evm and signal acquistion
specs.evm=mxx;
specs.bn=bnxx;
spec.ta=ta0;
specs.ua=ua0;
specs.mf = using_mf;

specs.xoc=xocc;
specs.cr=cr;
specs.crvec.pre=crvec(2)+38400;
specs.crvec.all=crvec(1)+38400;
specs.crvec.variation.pre=dv2(2);
specs.crvec.variation.all=dv(1);
specs.fs=fs;
specs.t0=t0;
specs.bits=bits;% payload bits
specs.pn.i=pni; % i chips is hex
specs.pn.q=pnq; % q chips in hex
specs.pn.posi = posi; % chip error positions
specs.pn.posq = posq; % chip error positions
specs.evmo=evmo;
if(exist('evmbr'))
    specs.evmbr = evmbr;
    specs.snrbr=snrbr;

end
specs.evmb=evmb;
specs.snrb=snrb;
specs.evmb_ss = evmb_ss;

specs.pkamp=pkamp;
specs.cno=cno;
specs.pn.invi=invi;
specs.pn.invq=invq;
specs.pn.ierr=fix(fchip-pki);
specs.pn.qerr=fix(fchip-pkq);
specs.mode=mode;
specs.fchip=fchip;
specs.p2p = p2p;
specs.rft.rt = rt;
specs.rft.ft = ft;
specs.rft.dur = dur;
