% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: qpskdemod 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: qpskdemod.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

%
%
%   The function doit will load a burst of samples 
%   and run through the T21 measurements
% 
function z=qpskdemod(y)

ini_phase=0;

% create the modmap
modmap = zeros(1,4);
modmap(1) = cos(ini_phase + pi/4) + 1i* sin(ini_phase + pi/4);
modmap(2) = cos(ini_phase + 7*pi/4) + 1i * sin(ini_phase + 7*pi/4);
modmap(3) = cos(ini_phase + 3*pi/4) + 1i * sin(ini_phase + 3*pi/4);
modmap(4) = cos(ini_phase + 5*pi/4) + 1i * sin(ini_phase + 5*pi/4);

z = genqamdemod(y,modmap);
