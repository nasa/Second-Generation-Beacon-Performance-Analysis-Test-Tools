% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: bin2hex 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: bin2hex.m 11 2019-09-23 13:10:04Z reesebo $
% ///            %
%
% function bin2hex( message )
%   message = binary bit vector, must not be string
%
% also see hex2bin and str2bin
%

function o=bin2hex(msg,varargin)

if(nargin>1)
    front=varargin{1};
else
    front=0;
end
sz=size(msg);

ln=length(msg);

if(ln<4)
    o=dec2hex(bin2dec(num2str(msg)));
    return;
end

for(ix=1:sz)
    message=msg(ix,:);
    begin=0;
    m = mod(length(message),4);
    endcap=[];
    if(m>0)
        if(~front)
            x=message(1:m);
            message = message(m+1:end);
            s=sprintf('%d',x);
            oo(1)=dec2hex(bin2dec(s));
            begin=1;
        else          
            message=[message zeros(1,4-m)];
            begin=0;
        end
    end
    y=reshape(message,4,length(message)/4);

    for(i=1:length(message)/4)
        x=y(:,i);
        s=sprintf('%d%d%d%d',x(1),x(2),x(3),x(4));
        oo(begin+i)=dec2hex(bin2dec(s));
    end
    o(ix,:) = [oo endcap];
end