% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: halfsinepulse 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: halfsinepulse.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

function [matched] = halfsinepulse(x,fsi,fchip)

sampleRate = (fsi/fchip);
hsp = sin(0:pi/sampleRate:(sampleRate)*pi/sampleRate);
decimationFactor = 1;
matchedFilter = dsp.FIRDecimator(decimationFactor,hsp);
matched = matchedFilter(x);
end

