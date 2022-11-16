% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: getsgbid
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: getsgbid.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
% 		MSB fixed bit:                              1
%       10-bit binary for 367 USA Country Code:		010 1101 111
% 		3 Fixed bits:                               1 01
% 		16-bit TAC:                                 0000 0000 0000 0001
% 		14-bit Beacon Serial Number:                00 0001 0000 0001
%       1-bit Test Protocol Flag                    0
% 		3-bit Aircraft/Vessel ID Type:              100
% 		44-bit Aircraft/Vessel ID:                  0 0000 0011 1001 1101 1011 0000
%                                                   0000 0000 0000 0000 000
%msg='00001004E3452FC7C043781D7559771A3A7FFFC56C7A300D1A835F7310870F80';
function o=getsgbid(hex)

bits=str2bin(hex2bin(hex,250,1));

tac=bits(1:16);
sn=bits(17:30);
cc=bits(31:40);
test=bits(43);
vidtype=bits(91:93);
vid=bits(94:137);

bid=[1 cc 1 0 1 tac sn test vidtype vid];

o=bin2hex(bid);