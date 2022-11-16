% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: doit 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 16 $ $Date: 2022-10-14 10:17:46 -0400 (Fri, 14 Oct 2022) $
% ///				$Id: doit.m 16 2022-10-14 14:17:46Z reesebo $
% ///            

%
%
%   The function doit will load a burst of samples 
%   and run through the T21 measurements
%   

%doit
clear;

%-----------User Inputs
pcal=53;  % power calibration setting.  If an attenuator/coupler is used, the value of the attenuation is input here.

% enter a date/time to filter the set of input waveform files (optional)
%dx=datenum('2/21/2021');
%enter the directory of matlab burst files from RSA.
xdir ='Test\';
% output matlab filename to save specs for post analysis if desired
outfn=['usernamehere' datestr(now,'yyyy_mm_dd_hh.MM.ss') '.mat'];


%-----------End User Inputs


f=dir([xdir '\*.mat']);
if(exist('dx') && ~isempty(dx))
    tx = [f.datenum];
    f=f(find(tx>dx));
end
global outdir;

for(ix=1:length(f))
    fprintf('loading: %s...', f(ix).name);
    load([xdir '\' f(ix).name]); 
    fsi=1/XDelta; % rsa sample time
    x=double(Y);  % rsa saves in single, convert to double
    fprintf('complete\n');
    
    risetime=0.04;  % padding for beginning of message.
    fchip=38400;    % ideal chip rate

    clear specs;
    
    % prepare file names and directories
    mkdir([xdir '\' f(ix).name(1:end-4)]);
    outdir=[xdir '\' f(ix).name(1:end-4)];
    csvfile = [outdir '\' f(ix).name(1:end-4) '.csv'];
    xlsxfile = [outdir '\' f(ix).name(1:end-4) '.xlsx'];
    outfile = [outdir '\' outfn];
    %close all;

    % processing steps:
    [specs] =getpayload(x,fsi,risetime,pcal);

    specs.filename = f(ix).name;

    if(specs.error == false)
        [specs]=signalquality(x,fsi,specs);
    end
    
    if(specs.error == false)
        [specs]=offsetcheck(x,fsi,specs);
    end
    
    if(specs.error == false)
        [specs]=testpayload(specs);
    end

    if(specs.error == false)
            specs=analyzemessage(specs);
    end
    %close all;

    
    if(specs.error == false)
        allspecs(ix) = specs;
        
        tab=specs.ms.show;

        fprintf('write measurement results to %s...', xlsxfile);
        writemeasurements(specs,xlsxfile);
        fprintf('complete\n');
        
        fprintf('write payload results...');
        writetable(tab,csvfile,'writerownames',true,'FileType','text','Delimiter',',','QuoteStrings',true);
        fprintf('complete\n');


    else
        fprintf('error message: %s\n', specs.errormessage);
    end
end
% save all specs for post analysis
save(outfile,'allspecs');
disp('Done!');