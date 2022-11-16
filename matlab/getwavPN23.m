% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: getwavPN23
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: getwavPN23.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
function [obm, xs]=getwavPN23(filename, mode)

if(~exist('mode','var'))
    mode=0;
end
fchip=38400;
sampsPerSym=4;

[ gci, gcq ] = getPN23( fchip, mode );

% % set data
canned={ '01234567', '89abcdef', 'deadbeef', '01234567', '89abcdef', 'deadbeef', 'FFFc2c15', '9136dcc0'};

arm=[];
for(ix=1:length(canned))
    arm=[arm dec2bin(hex2dec(canned{ix}),32)];
end

data=[str2bin(arm(1:250))];


[xs, obm] = oqpsk_spread(data, gci, gcq, filename);

% basis_matlab =makewaves_matlab(sampsPerSym);
% 
% obm = oqpsk_spread_basis(xint, basis_matlab,sampsPerSym);% using basis from arm 

%obm=obm/max(abs(obm));
