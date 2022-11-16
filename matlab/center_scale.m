% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: center_scale 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: center_scale.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

function [xb,cf,cr]=center_scale(xb,fchip,fs,fig)

    xb=xb(:).';
    shift=fix(fs/fchip/2);

    xshift=[xb(1+shift:end) zeros(1,shift)];
    
    xx=xshift.*xb;
if(fig)
    sfigure(fig);clf;
    pltstr='';
else
    pltstr=0;
end
    [~, cfr, ~, ~, ~]=specpeak_sr0(xx,-fs/2-1,fs/2+1,fs,pltstr); % (1:fix(167e-3*fs)) preamble cr
    cf1=max(cfr(1:2,1));
    cf2=min(cfr(1:2,1));

    cr=(cf1-cf2)/2;

    cf=(cf1-cr)/2;

    grid;
    xb=xb.*exp(-2*pi*1j*cf/fs*[0:length(xb)-1]);
    %xb=xb/max(abs(xb))/0.7071;

    xshift=[xb(1+shift:end) zeros(1,shift)];
    xx=xshift.*xb;
    
if(fig)
    pltstr=sprintf('cf=%f,cr=%f', cf, cr);
else
    pltstr=0;
end
    specpeak_sr0(xx,0,fs/2,fs,pltstr);