% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: sgb_bch_encode.m 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: sgb_bch_encode.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

function y = sgb_bch_encode( m )
% m is the binary message vector
% note n,k are not used here...the generator polynomial is created
% elsewhere using g=bchgenpoly(n,k);
% g is the generator polynomial
g=[1 1 1 0 0 0 1 1 1 1 1 1 0 1 0 1 1 1 0 0 0 0 1 0 1 1 1 0 1 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 1 0 1 1 1];

% get vector lengths
LM = length(m);
LG = length(g);
LY = LM + LG - 1;

% initialize vectors
y = zeros(1,LY);
y(1:LM) = m;
x = y;

% find the remainder
idx = 1:LG;
for k=1:LM
    if ( x(k) == 1 )
      x(idx) = xor( x(idx), g );
    end;
    idx = idx + 1;
end;

% attach remainder
y(LM+1:end) = x(LM+1:end);