% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: specpeak_sr2 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: specpowr_sr2.m 11 2019-09-23 13:10:04Z reesebo $
% ///            


function [centerfreq, cfr, CN0, CPow, N0]=specpowr_sr2(x,min_f,max_f,fs,pcal,str)

cfr=0;

if(str==0)
    plotstuff=0;
else
    plotstuff=1;
end
polypeak=1;

n = length(x);
N = 2.^(fix(log10(n)/log10(2)));
N = min(N,2^17);
faxis = [-N/2:N/2-1]*fs/N;

if(1)

%since we are not using a window function, the peak is sharper, but we do
%have spectral leakage that we are not accounting for.
% todo: add window and 3 point peak interpolation


x=x(1:N).^4; % 4th power for spread signals

W = blackman(N);
if(size(W,1) ~= size(x,1))
    W=W';
end
X = fft(x.*W/sum(W),N);

X = abs(fftshift(X));

else
    ino = fix(length(x)/Ni);
    for(i=1:ino)
        Xm(i,:)=fft(x(Ni*(i-1)+1:Ni*i),N);
    end
    X=fftshift(abs(mean(Xm)));
end

minb=find(faxis>min_f);
maxb=find(faxis<max_f);
minb=minb(1);
maxb=maxb(end);

try
    tryit=100;
    while(tryit)
        [maxtab,mintab]=peakdet(X(minb:maxb),tryit/100*max(X(minb:maxb)));
        if(size(maxtab,1) > 1)
            break;
        end
        tryit=tryit-10;
    end
    maxtab=flipud(sortrows(maxtab,2));
    cfr=maxtab;
    cfr(:,1) = faxis(cfr(:,1)+minb-1);
    pk=maxtab(1,2);

    centerfreq=faxis(maxtab(1,1)+minb-1);

catch
    centerfreq=faxis(minb-1);
end

if(polypeak)% interpolate peak here...
    try
        cfx = maxtab(1,1)+minb-1;
    catch
        cfx=minb-1;
    end
    %px=faxis(cfx-1:cfx+1)';
    px=faxis(N/2:N/2+2)';
    py=X(cfx-1:cfx+1);
    px=px(:);
    py=py(:);
    coeff=polyfit(px,py,2);
    pki=-coeff(2)/(2*coeff(1));
    centerfreq=pki+faxis(cfx);
    pk=polyval(coeff,pki);

    
else
[pk,pki]=max(abs(X(minb:maxb)));
%centerfreq = fs/N * (minb+pki);
centerfreq = fs/N * (minb+pki-1);
end



%cfr=peak
noise = 10*log10(median(X));
N0 = noise - 10*log10(fs/N);
CN0=10*log10(pk) - (noise - 10*log10(fs/N));

CPow = 10*log10(pk)/2; % 4th power needs to be reduced by 2
CPow = CPow-10*log10(100)+30+0.5+pcal;% convert to power, convert to dBm


if(plotstuff)
    %sfigure(2);
    plot(faxis,10*log10(X));
    hold on;
    plot((centerfreq)*[1 1],[noise 20*log10(abs(pk))],'y')
    plot(centerfreq,10*log10(abs(pk)),'go');
    %plot(cfr(:,1),20*log10(cfr(:,2)),'go');
    plot([faxis(1) faxis(end)], noise*[1 1],'r');
    hold off;
    title(sprintf('Power = %0.2fdBm @ %f %s\n',CPow, centerfreq, str));
    xlabel('Frequency(Hz)');
    ylabel('Amplitude(dB)');
    %axis([min_f max_f sort([noise CPow])]);
    drawnow;
end