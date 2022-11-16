% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: hex2bin 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: hex2bin.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
function o=hex2bin(x,numbits,varargin)

        if(nargin>2)
            front=varargin{1};
        else
            front=0;
        end

       
        o=[];

        for i=1:length(x)
            temp = sprintf('f%s', x(i));
            temp = dec2bin(hex2dec(temp));
            o = [ o temp(5:end)];
        end

        if(nargin>3)
            rightshift=varargin{2};
            if(rightshift>0)
                o=o(1:end-rightshift);
            elseif(rightshift<0)                  
                o=[o(abs(rightshift)+1:end) 0 0];
            end
        end
        
        if(~front)
            outdelta = length(o) - numbits+1;
            if(outdelta > 0)
                o = o(outdelta:end);
            end
        else
            if(length(o) >numbits)
                o=o(1:numbits);
            else
                o=[num2str(zeros(1,numbits-length(o))) o];
                o=o(o~=' ');
            end
        end
     end