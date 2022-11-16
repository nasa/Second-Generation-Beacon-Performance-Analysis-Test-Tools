% ///  Copyright(c) 2017 by Reese Bovard, all rights reserved
% ///  Proprietary Information of Concentric Real Time, LLC
% ///  Subject to Non-Disclosure.  Do Not Copy or Distribute
% ///  
% ///     Project: Second Generation Software Defined Receiver for Cospas-Sarsat DSSS Beacons
% ///     
% ///     Module: docorrb 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision$ $Date$
% ///				$Id$
% ///            
function [tdoa, pk, inv] = docorrb(fx,fy,fs,plotstuff)

if(plotstuff)
    sfigure(plotstuff);clf;
end
tx=1:length(fx);
ty=1:length(fy);
tcom=0;

if(length(tx) > length(ty))
    t=tx;
    
else
    t=ty;
end

w=t(end);
pxh = fx;
pyh=fy;


if(plotstuff)
% fprintf('x=%20.20f vs y=%20.20f)\n', tcom-tx(1), tcom-ty(1));

%figure(88);
subplot(311);plot(tx-tcom,pxh);
title('signal');
ylabel('amplitude');
xlabel('time');
%axis([0 w -2000 2000]);
%hold on;plot(ones(1,length(t)*t-tcom, -2000:2000,'r');hold off;
%title(sprintf('x=%s',fx));

subplot(312);plot(ty-tcom,pyh);
hold on;plot([0 0],[ -4 4],'r');hold off;
title('sync reference and signal');
ylabel('amplitude');
xlabel('time');
%axis([0 w -2000 2000]);
%hold on;plot(ones(1,length(-2000:2000))*t-tcom, -2000:2000,'r');hold off;
%title(sprintf('y=%s',fy));
%pause;
end

% 3. xcorr complex
[c,lags] = xcorr(pxh,pyh);
% 4. abs of output is lag: associated 'x' timetag (as long as we are
% consistent)
C=abs(c);

[kpeak,kpeaki]=max(C);
%     px=lags(kpeaki-1:kpeaki+1);
%     py=C(kpeaki-1:kpeaki+1);py=py/kpeak;
%     coeff=polyfit(px(:)',py(:)',2);
%     pki=-coeff(2)/(2*coeff(1));
%     pk=polyval(coeff,pki)*kpeak;
%     %plot(lags,C);
%     rng = pki-3:pki+3;
%     pkrng=polyval(coeff,rng)*kpeak;
pk=kpeak;
pki=lags(kpeaki);

offset = lags(kpeaki); % zero based
[klow,klowi]=min(c);

if(kpeaki==klowi && abs(klow)==kpeak)
    inv=true;
else
    inv=false;
end
%offset=length(pxh)-mcmi;
%offset=mcmi-length(t);
tdoa=round(offset);
%t=tcom;

if(plotstuff)
subplot(311);hold on;plot([tdoa tdoa],[ -4 4],'r');hold off;
ylabel('amplitude');
%subplot(312);hold on;plot([tdoa tdoa],[ -4 4],'r');hold off;
ylabel('amplitude');
xlabel('time');

%subplot(313);plot(lags,10*log10(C));hold;plot(mcmi,10*log10(pk),'ro');hold off;
subplot(313);plot(lags,(C));hold;plot(pki,(pk),'ro');hold off;
ylabel('amplitude');
xlabel('lags');

%title(sprintf('x has a %f vs ix=%f sec lag relative to %f (x=%f,y=%f)',tdoa, (x-y), t, tx(1), ty(1)));
title(sprintf('x has a %d sample (%02.04g sec) lag relative to y pk=%d',tdoa,tdoa/fs,pk ));
%pause;
%drawnow;
%figure(89)
% if(tdoa<0)
%     tdoa=0;
% end
subplot(312);
if(tdoa>0)
plot((pxh(fix(tdoa)+1:end))/max(abs(pxh(1:end))),'b'); hold on;
plot((pyh(1:end))/max(abs(pxh(1:end))),'r');
elseif(tdoa<0)
    plot((pxh(1:end))/max(abs(pxh(1:end))),'b'); hold on;
    plot((pyh(abs(tdoa)+1:end))/max(abs(pxh(1:end))),'r');
else
    plot((pxh(1:end))/max(abs(pxh(1:end))),'b'); hold on;
    plot((pyh(1:end))/max(abs(pxh(1:end))),'r');
end
plot([tdoa tdoa],[ -4 4],'r');
hold on;
% plot(imag(pyh)/max(abs(pyh)),'r');
% plot(real(pyh)/max(abs(pyh)),'r');

hold off;
%figure(88);
end

% 5. deltaphase/dt of adjacent sample from peak will be freq
%f=(angle(c(mcmi))-angle(c(mcmi-1)))* fs;
f=0;

% fprintf('offset=%f, tdoa=%f, tref=%f, freq=%f\n', offset, tdoa, t, f);