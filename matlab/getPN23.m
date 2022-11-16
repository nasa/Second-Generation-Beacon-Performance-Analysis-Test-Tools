% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: getPN23 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: getPN23.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
% Copyright © 2014 Concentric Real Time, LLC. 
% CRT reserves all rights in the Program as delivered.  The Program or 
% any portion thereof may not be reproduced in any form whatsoever except 
% as provided by license.  
% This software was developed for NASA under Contract NNG13CR48C.
function [gci,gcq]=getPN23(N,mode)

if(~isempty(mode) && strcmp(mode,'selftest'))
    selftest=2;
else
    selftest=0;
end

init={...
... %MSB                        LSB
    '0000 0000 0000 0000 0000 001';...
    '0011 0101 1000 0011 1111 100'  ;
    '1010 0101 1001 0011 1110 000'; ...
    '0111 1001 1101 0010 0101 000'; ...
    };
     % x^23      x^18                                x^0
taps= [1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
win=selftest+(1:2);
for(ix=win)
    A=init{ix,:};
    initb = str2bin(A(find(~isspace(A))));
    hpn0 = comm.PNSequence('Polynomial', taps, ... % [23 18 0]
          'SamplesPerFrame', N, 'InitialConditions',initb);
    x1 = step(hpn0);

    if(mod(ix,2)==1)
        gci=x1;
    else
        gcq=x1;
    end
     tapv=taps(2:end);  % remove x^23
%     x2 = mls_mat(tapv, initb,N);

    %disp([ A ': ' bin2hex(x1(1:N)) ':' bin2hex(x2)]);
    %disp([ A ': ' bin2hex(x1(1:N))]);

    
end