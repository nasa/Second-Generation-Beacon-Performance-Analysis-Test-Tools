% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: signalquality 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: signalquality.m 15 2022-09-29 15:45:13Z reesebo $
% ///            

function [specs]=signalquality(x,fsi,specs)

if(strcmp(specs.mf,'sin') == 1)
    specs.sq2.offset=[];
    specs.sq2.dev=[];
    specs.sq2.dev100=[];
    specs.sq2.dev167=[];
    specs.sq2.s2=[];
    specs.sq2.fs=[];
    specs.sq2.maxMinusMin = [];
    return
end

t0=specs.t0;
fchip=specs.cr;
bn=specs.bn;

fchip=fix(fchip);

sps=40;
fs=fchip*sps;

fprintf('resample...');
[a,b]=rat(fs/fsi);
xo=resample(x,a,b);
xo=xo/max(abs(xo));
fprintf('complete\n');
fig=200;


fprintf('profile power and frequency over burst...');

t0s=fix(t0*fs);
fchip=38400;
fs=sps*fchip;
if(t0s<0)
    t0s=1;
end

span=t0s:t0s+fix(fs)-1;

if(span(end)> length(xo))
    span=1:length(xo);
end
xob=xo(span); 

[xoc,offset,~]=center_scale(xob,fchip,fs,0);

fprintf('complete\n');
b = [ 1 2 3 4 3 2 1 ]; b = b/sum(b);
y0 = filter( b, 1, xoc );

fprintf('correct carrier phase...');
[y0, p0, t0]=carrierphasecorrection(y0,fs,bn,fchip,fig);
%[y0 ,p0, t0]=carrierphasecorrection2(xoc,fs,sps,2,fig);p0=p0(:).';
saveasrsb(fig,'carrierphasecorrect.png');fig=fig+3;
[dev_100,s2,maxMinusMin]=shorttermcarrier_fgblike(p0,t0,fig);fig=fig+1;
[dev_167,o2,maxMinusMin]=shorttermcarrier(p0,t0,167.7e-3,fig);
[dev_42,o2,maxMinusMin]=shorttermcarrier(p0,t0,41.66667e-3,fig);
saveasrsb(fig,'shorttermcarrier.png');fig=fig+1;
fprintf('complete\n');

specs.sq2.offset=offset+o2;
specs.sq2.dev=dev_42;
specs.sq2.dev100=dev_100;
specs.sq2.dev167=dev_167;
specs.sq2.s2=s2;
specs.sq2.fs=fs;
specs.sq2.maxMinusMin = maxMinusMin;
