% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: dospread 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: dospread.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
function [ sii, sqq ]=dospread(bits,mode)
% mode='selftest' or 'normal';
% build waveform NRZ with lookup table using sinc pulse

fchip=38400;
PL = 6400;
N  = 300;

[ si, sq ] = getPN23( fchip, mode );


% apply modulation to sequence
siq = zeros(1,2*fchip+2);
siq(1:2:end-3) = si;
siq(2:2:end-2) = sq;
siq(end-1) = si(end);
siq(end)   = sq(end);

miq = siq;
r = [ 1 : 2 :512 ];
for n=0:N/2-1
    miq(r  ) = xor( miq(r  ), bits(2*n+1) ); 
    miq(r+1) = xor( miq(r+1), bits(2*n+2) ); 
    r = r + 512;
end;

% modulated and spread
sii=miq(1:2:end-2);  
sqq=miq(2:2:end-1);