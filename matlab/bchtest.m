% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: bchtest 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: bchtest.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

function [outbits,err,errpos]=bchtest(bits)

bits=bits(:);
N=250;
K=202;
dec = comm.BCHDecoder('CodewordLength',N,'MessageLength',K,'NumCorrectedErrorsOutputPort',true);

[outbits,err]=step(dec,bits);
outbits=outbits(:);
if(err>0)
    fprintf('BitErrors Detected\n');
end
errpos = find(xor(bits(1:202), outbits) ==1);
fprintf('complete\n');
