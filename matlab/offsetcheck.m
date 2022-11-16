% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: offsetcheck 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: offsetcheck.m 15 2022-09-29 15:45:13Z reesebo $
% ///            

function specs=offsetcheck(x,fsi,specs)

t0=specs.t0;
fchip=specs.cr;
mode=specs.mode;
data=specs.bits(51:end);
bn=specs.bn;

fig=300;


fchip=round(fchip);

sps=100;
fs=fchip*sps;

[a,b]=rat(fs/fsi);
xo=resample(x,a,b);
xo=xo/max(abs(xo));
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
fprintf('recenter frequency...');

[xoc,~,~]=center_scale(xob,fchip,fs,0);

fprintf('complete\n');

if(strcmp(specs.mf,'sin') ==1 )
    y0=xoc;
    bn=0.01;
    y0=carrierphasecorrection2(y0,fs,sps,bn,3,4);
else
    b = [ 1 2 3 4 3 2 1 ]; b = b/sum(b);
    y0 = filter( b, 1, xoc );

    fprintf('correct carrier phase...');
    [y0]=carrierphasecorrection(y0,fs,bn,fchip,0);
end

[xm, ~]=getwfPN23_upsamplewithdata(sps,'',specs.mf, mode,data); 
fprintf('complete\n');

% check realtive shift
fprintf('check offset...');

sfigure(fig);clf; hold on;

xr=real(xm);
xi=imag(xm);

[cm0, nxr]=xcorr(xr,y0);
[cm1, nxi]=xcorr(xi,y0);

cm0=abs(cm0)/max(abs(cm0));
cm1=abs(cm1)/max(abs(cm1));

pki=polypeak(cm0,nxr);
pkq=polypeak(cm1,nxi);

chiperror = abs((pki-pkq)/fs*fchip*100);

specs.iqoffset=chiperror;

title(sprintf('I vs Q Offset Error = %0.2f%% ',chiperror));
legend('I Correlation Fit','I Peak','Q Correlation Fit','Q Peak')
xlabel('offset error (samples)');
ylabel('Correlation Amplitude');
saveasrsb(fig,'iqoffset.png');
fprintf('complete\n');

% find quadratic interpolated peak
function pk=polypeak(cm,nx)

[~,ix]=max(cm);
yf=cm(ix-1: ix+1);
xf=nx(ix-1: ix+1);
xfp=[xf(1):.01:xf(end)];
p=polyfit(xf,yf,2);
pv=polyval(p,xfp);
plot(xfp,pv);
pk=-p(2)/(2*p(1));
plot(pk,polyval(p,pk),'x');