% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: docno 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: docno.m 15 2022-09-29 15:45:13Z reesebo $
% ///
function [cno, evm, snr]=docno(s)

ss=[real(s) imag(s)];

o=ss(ss>0);
z=ss(ss<0);

r = (mean(o.^2) + mean(z.^2));
v = var(o) + var(z);

evm=sqrt(v/r);
snr=-20*log10(evm);
cno=snr+10*log10(38400*2);
