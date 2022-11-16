% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: t21_test 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: t21_test.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
% ///
% ///
% ///   The function t21_test will load each burst from specified directory  
% ///   xdir and run through the T21 measurements
% ///   

function t21_test(xdir,pcal)

f=dir([xdir '\*.mat']);

for(ix=1:length(f))
    load([xdir '\' f(ix).name]); 
    fprintf('Output: %s\n', [xdir '\' f(ix).name]);
    risetime=0.04; % seconds
    fsi=1/XDelta;
    x=double(Y);
    fchip=38400;

    mkdir([xdir '\' f(ix).name(1:end-4)]);
    outdir=[xdir '\' f(ix).name(1:end-4)];
    csvfile = [xdir '\' f(ix).name(1:end-4) '.csv'];
    xlsxfile = [xdir '\' f(ix).name(1:end-4) '.xlsx'];
    close all;

    [specs] =getpayload(x,fsi,risetime,pcal);

    [specs]=signalquality(x,fsi,specs);

    [specs]=offsetcheck(x,fsi,specs);
    
    [specs]=testpayload(specs);

    specs=analyzemessage(specs);

    tab=specs.ms.show;

    compilespecs(tab, specs, xlsxfile);

    writetable(specs.ms.show,csvfile,'writerownames',true,'FileType','text','Delimiter',',','QuoteStrings',true);

end

disp('Done!');