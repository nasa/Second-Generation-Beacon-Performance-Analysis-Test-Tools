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
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: getdetect.m 15 2022-09-29 15:45:13Z reesebo $
% ///            
function [xocc, t0, pkamp]=getdetect(xoc,fs,fchip,mode)
fig=700;

pkamp=0;

if(1)
fprintf('detect signal...');
maxshift=fix(200e-3*fs);
sfigure(fig);fig=fig+1;
crossings=detectpn23(xoc,maxshift,mode,fs/fchip,fig);

if(isempty(crossings) || (abs(crossings(1,3)) > 1000 && isempty(abs(crossings(:,3)) <1000)) )
    xoc=conj(xoc);
    crossings=detectpn23(xoc,maxshift,mode,fs/fchip);
    if(~isempty(crossings))
            title('Detection, Inverted Spectrum');
    else
        disp('Error! Signal Not Detected!!');
        t0=0;
        xocc=[];
        return;
    end
else
    title('Detection, Non-inverted Spectrum');
end
[pkamp,ix]=max(crossings(:,2));
Nd=crossings(ix,1);
t0 = Nd/fs;
fprintf('complete\n');
len=min([length(xoc(Nd:end))-1 fix(fs*1.01)]);
xocc=xoc(Nd:Nd+len);

else
[off0,pk0,inv0]=docorrb(real(xoc),real(xs),fchip,1);
[off1,pk1,inv1]=docorrb(real(xoc),imag(xs),fchip,1);

if(pk1>pk0)
    xoc=conj(xoc(off1:end));
    t0=t0+off1/fs;
else
    xoc=(xoc(off0:end));
    t0=t0+off0/fs;
end
end