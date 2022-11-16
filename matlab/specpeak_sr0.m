% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: specpeak_sr0 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: specpeak_sr0.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

function [centerfreq, cfr, CN0, CPow, N0]=specpeak_sr0(x,min_f,max_f,fs,str)

cfr=0;

if(str==0)
    plotstuff=0;
else
    plotstuff=1;
end
polypeak=1;

%N=64*1024*2;
N=length(x);
N=N-mod(N,2);
Nl=N;
% if(Nl<N)
%     N=Nl;
% end

faxis = [-N/2:N/2-1]*fs/N;

if(1)

%since we are not using a window function, the peak is sharper, but we do
%have spectral leakage that we are not accounting for.
% todo: add window and 3 point peak interpolation
W=blackman(N);
X = fftshift(abs(fft(x(1:N).*W',N)));

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
    for(ix=1:length(maxtab(:,1)))
        
        try
            cfx = maxtab(ix,1)+minb-1;
        catch
            cfx=minb-1;
        end
        px=faxis(N/2:N/2+2)';
        py=X(cfx-1:cfx+1);
        px=px(:);
        py=py(:);
        coeff=polyfit(px,py,2);
        pki=-coeff(2)/(2*coeff(1));
        cfr(ix,1) =pki+faxis(cfx);
        pk=polyval(coeff,pki);
        cfr(ix,2) = pk;
    end
    centerfreq=cfr(1,1);
    pk=cfr(1,2);
else
[pk,pki]=max(abs(X(minb:maxb)));
centerfreq = fs/N * (minb+pki-1);
end



%cfr=peak
noise = 20*log10(median(X));
N0 = noise - 10*log10(fs/N);
CPow = 20*log10(pk);

CN0=20*log10(pk) - (noise - 10*log10(fs/N))+10*log10(N/Nl);

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
    %axis([min_f max_f sort([noise CPow])]);
    drawnow;
end