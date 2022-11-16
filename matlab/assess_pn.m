% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: assess_pn 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: assess_pn.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

function [isuccess,qsuccess,invi,invq,pki,pkq,ipositions,qpositions]=assess_pn(bits,pnbin,fchip,mode)
fig=400;
% assess pn code:
inpn=zeros(1,2*fchip);
inpn(1:length(pnbin)) = pnbin;
rci=inpn(1:2:end);
rcq=inpn(2:2:end);

% need to get the modulated message with known data bits here.  Only the
% preamble matches up unless we have the known message.
[ si, sq ]=dospread(bits,mode);

[off0,pkin,invin]=docorrb(2*si(1:fchip)-1,2*rci-1,fchip,0);
[off1,pkii,invii]=docorrb(2*si(1:fchip)-1,2*rcq-1,fchip,0);
frntpad=[];
backpad=[];
if(pkin>pkii)
    pki=pkin;
    invi=invin;
    if(off0>0)
        frntpad=(zeros(1,off0));
    elseif (off0<0)
        off0=abs(off0);
        backpad=(zeros(1,off0));
    end
    pni = [rci(off0+1:end) backpad];
    pni=pni(1:fchip);
    if(invi)
        pni=~pni;
        rci=~rci;
    end
    [~,pki,invi]=docorrb(2*si(1:fchip)-1,2*pni-1,fchip,fig);fig=fig+1;
    misses = xor(si,pni);
    [~,ipositions]=find(misses>0);
    pki=38400-sum(misses);
    pni=bin2hex(pni,1);
    fprintf('%s\n',bin2hex(si,1));
    fprintf('%s\n',pni);
else
    pki=pkii;
    invi=invii;
    if(off1>0)
        frntpad=(zeros(1,off1));
    elseif (off1<0)
        off1=abs(off1);
        backpad=(zeros(1,off1));
    end
    pni = [rcq(off1+1:end) backpad];
    pni=pni(1:fchip);
    if(invi)
        pni=~pni;
        rcq=~rcq;
    end
    [~,pki,invi]=docorrb(2*si(1:fchip)-1,2*pni-1,fchip,fig);fig=fig+1;
    misses = xor(si,pni);
    [~,ipositions]=find(misses>0);
    pki=38400-sum(misses);

    pni=bin2hex(pni,1);
    fprintf('%s\n',bin2hex(si,1));
    fprintf('%s\n',pni);
    disp('swap');
end
[off0,pkqn,invqn]=docorrb(2*sq(1:fchip)-1,2*rcq-1,fchip,0);
[off1,pkqi,invqi]=docorrb(2*sq(1:fchip)-1,2*rci-1,fchip,0);
if(pkqn>pkqi)
    pkq=pkqn;
    invq=invqn;
    if(off0>0)
        frntpad=(zeros(1,off0));
    elseif (off0<0)
        off0=abs(off0);
        backpad=(zeros(1,off0));
    end
    pnq = [rcq(off0+1:end) backpad];
    pnq=pnq(1:fchip);
    if(invq)
        pnq=~pnq;
        rcq=~rcq;
    end
    [~,pkq,invq]=docorrb(2*sq(1:fchip)-1,2*pnq-1,fchip,fig);fig=fig+1;
    misses = xor(sq,pnq);
    [~,qpositions]=find(misses>0);
    pkq=38400-sum(misses);
    pnq=bin2hex(pnq,1);
    fprintf('%s\n',bin2hex(sq,1));
    fprintf('%s\n',pnq);
else
    pkq=pkqi;
    invq=invqi;
    if(off1>0)
        frntpad=(zeros(1,off1));
    elseif (off1<0)
        backpad=(zeros(1,abs(off1)));
        off1=abs(off1);
    end
    pnq = [rci(off1+1:end) backpad];
    pnq=pnq(1:fchip);
    if(invq)
        pnq=~pnq;
        rci=~rci;
    end
    [~,pkq,invq]=docorrb(2*sq(1:fchip)-1,2*pnq-1,fchip,fig);fig=fig+1;
    misses = xor(sq,pnq);
    [~,qpositions]=find(misses>0);
    pkq=38400-sum(misses);
    pnq=bin2hex(pnq,1);    
    fprintf('%s\n',bin2hex(sq,1));
    fprintf('%s\n',pnq);
    disp('swap');
end
isuccess=pni;
qsuccess=pnq;

fprintf('Total number of matching chips in I: %d (Inverted: %d)\n', fix(pki), invi);
fprintf('Total number of matching chips in Q: %d (Inverted: %d)\n', fix(pkq), invq);
