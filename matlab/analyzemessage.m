% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: analyzemessage 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: analyzemessage.m 15 2022-09-29 15:45:13Z reesebo $
% ///            

function specs=analyzemessage(specs)

bits=specs.bits;

%npayload=202;
fprintf('parse message...');

ms=messageStruct;
ms.setvecb(bits(51:end));
fprintf('complete\n');

specs.ms=ms;
