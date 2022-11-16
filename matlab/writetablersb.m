% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: writetablersb 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: writetablersb.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
function writetablersb(tab,xlsxfile)
rownames = tab.Properties.RowNames;

for(ix=1:length(rownames))
    xlswrite(xlsxfile, rownames(ix), 'Payload', sprintf('A%d', ix+1));
end

cols=tab.Properties.VariableNames;

for(ix=1:length(cols))
    xlswrite(xlsxfile, cols(ix), 'Payload', sprintf('%s1', ix+65));
end

sz=size(tab);

for(ix=1:sz(2))
        if(~isnumeric(tab{1,ix}))
            xlswrite(xlsxfile, tab{:,ix}, 'Payload', sprintf('%s%d', ix+65,2));
        else
            xlswrite(xlsxfile, tab{:,ix}, 'Payload', sprintf('%s%d', ix+65,2));
        end
end