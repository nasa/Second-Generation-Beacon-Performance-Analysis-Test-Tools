% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: testpayload 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: testpayload.m 15 2022-09-29 15:45:13Z reesebo $
% ///            

function [specs]=testpayload(specs)

bits=specs.bits(:);
fprintf('message payload check using BCH...');

N=250;
K=202;
dec = comm.BCHDecoder('CodewordLength',N,'MessageLength',K,'NumCorrectedErrorsOutputPort',true);

[outbits,err]=step(dec,bits(51:end));
outbits=outbits(:);
if(err>0)
    fprintf('BitErrors Detected\n');
elseif(err<0)
        specs.error=1;
end
errpos = find(xor(bits(50+(1:202)), outbits) ==1);
fprintf('complete\n');

specs.payload.outbits=outbits;
specs.payload.err=err;
specs.payload.errpos=errpos;
specs.payload.uncorrected = bin2hex( [ 0 0 bits(51:end)']);


