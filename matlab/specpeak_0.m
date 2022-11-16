% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: specpeak_0 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: specpeak_0.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

function [centerfreq,pk]=specpeak_0(x,min_f,max_f,fs,str)

if(isempty(str) || str(1)==0)
    plotstuff=0;
else
    plotstuff=1;
end


N=64*1024;
Ni=32*1024;
Nl=length(x);

if(Nl<Ni)
    Ni=Nl;
end

faxis = [-N/2:N/2-1]*fs/N;


x=x(1:Ni);

%since we are not using a window function, the peak is sharper, but we do
%have spectral leakage that we are not accounting for.
% todo: add window and 3 point peak interpolation
X=fftshift(abs(fft(x,N))); 

minb=find(faxis>min_f);
maxb=find(faxis<max_f);
minb=minb(1);
maxb=maxb(end);


[y,maxtab]=max(X(minb:maxb));




    cfx = maxtab+minb-1;

    %px=faxis(cfx-1:cfx+1)';
    px=faxis(N/2:N/2+2)';
    py=X(cfx-1:cfx+1);
    px=px(:);
    py=py(:);
    coeff=polyfit(px,py,2);
    pki=-coeff(2)/(2*coeff(1));
    centerfreq=pki+faxis(cfx);
    pk=polyval(coeff,pki);

    




%cfr=peak
noise = 20*log10(median(X));
N0 = noise - 10*log10(fs/N);
CPow = 20*log10(pk);

CN0=20*log10(pk) - (noise - 10*log10(fs/N))+10*log10(N/Ni);


if(plotstuff)
    %sfigure(2);
    plot(faxis,20*log10(X));
    hold on;
    plot((centerfreq)*[1 1],[noise 20*log10(abs(pk))],'y')
    plot(centerfreq,20*log10(abs(pk)),'go');
    %plot(cfr(:,1),20*log10(cfr(:,2)),'go');
    plot([faxis(1) faxis(end)], noise*[1 1],'r');
    hold off;
    title(sprintf('C/N0 = %f @ %f %s\n',CN0, centerfreq, str));
    xlabel('Frequency(Hz)');
    ylabel('Amplitude(dB)');
    %axis([min_f max_f noise CPow]);
    drawnow;
end