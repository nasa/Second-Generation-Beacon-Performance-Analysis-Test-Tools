% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: matchfilter 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 16 $ $Date: 2022-10-14 10:17:46 -0400 (Fri, 14 Oct 2022) $
% ///				$Id: matchfilter.m 16 2022-10-14 14:17:46Z reesebo $
% ///            
function o=matchfilter(x,fs,fchip)

sps=fix(fs/fchip);
Hd=mf2(fs,sps-mod(sps,2));

o=filter(Hd.Numerator,1,x);

%quickfft(o,fs);