% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: chiprate3 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: chiprate3.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

%
%   [cr,amp]=chiprate3(xb,fs,fchip,shift,bounds,drw)  
%
%   Inputs:
%       xb: complex waveform samples
%       fs: sample rate of complex samples
%       fchip: expected chip rate of waveform
%       shift: the number of samples to shift a copy of the burst.
%       bounds: find the peak within the "bounds" window
%       drw: if true, will plot results, else will not
%
%   Outputs:
%       cr: single chip rate estimate 
%       amp: amplitude of the chip rate peak.
%
%   The function chiprate3 will estimate the chip rate over a given
%   portion of the chip rate domain using the given samples.
%   

function [cr,amp]=chiprate3(xb,fs,fchip,shift,bounds,fig)   

    xb=xb(:).';
    % shifted input vector
    xshift=[xb(1+shift:end) zeros(1,shift)];
    
    % chiprate domain
    %[x1,faxis]=quickfft0(xshift.*xb,fs);

    if(fig==0)
        str=0;
    else
        str='';
        sfigure(fig);clf;
    end
    [centerfreq, cfr, CN0, CPow, N0]=specpeak_sr0(xshift.*xb,-fs/2,fs/2,fs,str);
 
    if(size(cfr,1) > 2)
        [Y,I] =sort(abs(cfr),1,'descend');
        cfr=cfr(I(1:2,1),:);
    end
    
    try
%     cf1 = cfr(abs(cfr(:,1)-fchip) < 1,1);
%     cf2 = cfr(abs(cfr(:,1)+fchip) < 1,1);
    cf1=max(cfr(1:2,1));
    cf2=min(cfr(1:2,1));

    cr=(cf1-cf2)/2;
    catch
        cr=0;
    end
    if(fig ~= 0)
        title(sprintf('chip rate = %f', cr));
    end
    
