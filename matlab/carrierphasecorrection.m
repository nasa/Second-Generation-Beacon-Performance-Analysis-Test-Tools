% ///  Copyright(c) 2017 by Reese Bovard, all rights reserved
% ///  Proprietary Information of Concentric Real Time, LLC
% ///  Subject to Non-Disclosure.  Do Not Copy or Distribute
% ///  
% ///     Project: T21 Software Tools for Cospas-Sarsat DSSS Beacons
% ///     
% ///     Module: carrierphasecorrection 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 12 $ $Date: 2020-01-21 08:10:30 -0500 (Tue, 21 Jan 2020) $
% ///				$Id: carrierphasecorrection.m 12 2020-01-21 13:10:30Z reesebo $
% ///            

function [y,lamsav,t]=carrierphasecorrection(x,fs,bn,fchip,fig)

sps=round(fs/fchip);

N=fix(fs*0.4);

y=[];
%bn=0.01;
kp=sps; % sample per symbol
ko=2; % qpsk
psilast=0;
lam=0;
lamlast=0;
enlast=0;
theta=bn/(0.707+1/(4*0.707));
d=1+2*0.707*theta+theta^2;
gi=4*(theta^2/d)/(kp*ko);
gp=4*0.707*(theta/d)/(kp*ko);
for(n=2:length(x))    
    y(n)=x(n)*exp(-2*pi*j*lam);
    en=sign(real(y(n)))*imag(y(n)) - sign(imag(y(n)))*real(y(n));
    psi=gi*en+psilast;
    lam = gp*enlast+psilast+lamlast;
    psilast=psi;
    lamlast=lam;
    enlast=en;
    ensav(n)=en;
    lamsav(n)=lam;
end

lamsav=lamsav*2*pi; % phase error in radians

if(fig>0)
    sfigure(fig);clf;
    subplot(311);
    t=showiq(y,fs);
    title('corrected signal');
    ylabel('I/Q');
    subplot(312);
    plot(t,ensav);title('error (rads)');
    ylabel('\phi error');
    subplot(313);
    len=min([length(lamsav) fs]);
    plot(t,-lamsav); title('Phase Error');

    title('phase correction');
    xlabel('time (samples)');
    ylabel('phase err (rads)');


%     h=scatterplot(x,1,0,'b');hold on;
%     h=scatterplot(y,1,0,'yx',h);


    %hi=eyediagram(y(N:N+75*sps).',2*sps,sps,0,'b');
end

