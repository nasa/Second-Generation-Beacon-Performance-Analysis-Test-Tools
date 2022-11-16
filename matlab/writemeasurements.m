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
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: writemeasurements.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

%
%
%   The function doit will load a burst of samples 
%   and run through the T21 measurements
%   
function writemeasurements(specs, fn)

copyfile('SGB - US Results Template.xlsx', fn);

fields = { ...
    'k8', ... % rise time
    'k9' ,... % power before
    'k10',... % power after
    'k11',... % transmission time
    'k14',... % freq stab long
    'k16',... % freq stab short
    'k18',... % IQ errors
    'k20',... % chip rate pre
    'k22',... % chip rate pre vary
    'k24',... % chip rate all
    'k26',... % chip rate all vary
    'k28',... % IQ offset
    'k29',... % Pk to Pk amplitude
    'k30',... % EVM
    'k31',... % Spurious In
    'k32',... % Spurious Out
};

content= { ...
    specs.rft.rt*1000, ... % convert to ms
    0, ...
    0, ...
    specs.rft.dur*1000, ... % convert to ms
    specs.sq.offset/406.05, ... % convert to ppm
    specs.sq.dev, ... % in ppb
    specs.pn.ierr+specs.pn.qerr, ...
    specs.crvec.pre, ... % in chips/s
    specs.crvec.variation.pre, ... %in chips/sec^2
    specs.crvec.all, ... % in chips/s
    specs.crvec.variation.all, ... %in chips/sec^2
    specs.iqoffset/100, ... % in percent; excel expects ratio
    specs.p2p/100, ...  % in percent; excel expects ratio
    max(specs.evm)/100, ... % convert to ratio
    0, ...
    0 ...
    };

for(ix=1:length(fields))
    if(ischar(content{ix}))
        o=content(ix);
    else
        o=content{ix};
    end
    xlswrite(fn, o, 'Measure', fields{ix});
end



