% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: carrierphasecorrection2 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: carrierphasecorrection2.m 15 2022-09-29 15:45:13Z reesebo $
% ///            

function [y,frq,t]=carrierphasecorrection2(x,fs,sps,bn,type,debugit)


%persistent carrierSync;
lamsav=[];
ensav=[];
%sps=round(fs/fchip);
%     [mx]=max(abs(x));
%     x=x.*0.707/mx;

%N=fix(fs*0.4);

if(type==1)
    sfigure(102);clf;hold on;

    x=real(x)-mean(real(x)) + 1j* (imag(x)-mean(imag(x)));
    x1=x;
    idx=find(imag(x)<0);
    x1(idx) = -(x1(idx));
    plot(x1,'.');

    ph = mean(angle(x1));
    x1=x.*exp(-j*ph);
    plot(x1,'.');
    grid on; 
    drawnow;
    x=x1;
end

if(1)
    x=x(:);
    if(~exist('carrierSync') || isempty(carrierSync))
        switch(type)
            case 1
                type='BPSK';
            case 2
                type='OQPSK';
            case 3
                type='QPSK';
        end

        carrierSync = comm.CarrierSynchronizer( ...
        'SamplesPerSymbol',sps,'Modulation',type ,'NormalizedLoopBandwidth', bn, 'DampingFactor', .7071 );
    end
    %carrierSync.reset();
    [y,lamsav]=carrierSync(x);
else
    
y=[];
bn=0.01;
kp=sps; % sample per symbol
ko=2; % oqpsk,qpsk,bpsk=2
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
    if(type==2)
        en=sign(real(y(n)))*imag(y(n)) - sign(imag(y(n)))*real(y(n));
    else
        en=sign(real(y(n)))*imag(y(n));
    end
    psi=gi*en+psilast;
    lam = gp*enlast+psilast+lamlast;
    psilast=psi;
    lamlast=lam;
    enlast=en;
    ensav(n)=en;
    lamsav(n)=lam;
end


end
    len=min([length(lamsav) fs]);
    frq=-lamsav(1:len);
    t=(0:length(frq)-1)/fs;

if(debugit>1)
    fig=debugit;
    sfigure(fig);clf;
    subplot(311);
    showiq(y,fs);
    title('corrected signal');
    ylabel('I/Q');
    subplot(312);
    plot(ensav);title('error (rads)');
    ylabel('\phi error');
    subplot(313);
    plot(frq); title('freq');

    title('freq correction');
    xlabel('time (samples)');
    ylabel('frq err (Hz)');
end
% if(debugit>0)
%         h=scatterplot(x,1,0,'b');hold on;
%         h=scatterplot(y,1,0,'yx',h);
% 
%         %N=1;
%         %hi=eyediagram(y(N:N+75*sps).',2*sps,sps,0,'b');
%         hold on;
%         Constellation=.7071*[1+1i,1-1i,-1+1i,-1-1i];
%         SymbolMapping=[0,1,2,3];
%         % for jj=1:4
%         %         text(real(Constellation(jj))-0.15,...,
%         %         imag(Constellation(jj))+0.15,...
%         %         dec2base(SymbolMapping(jj),2,2));
%         %     plot(real(Constellation(jj)),...,
%         %         imag(Constellation(jj)),'rd');
%         % end
%         scatterplot(Constellation,1,0,'r+',h);
%         title('phase corrected symbol map');
%         axis([-1.5 1.5 -1.5 1.5]);
%         grid on;
%     
% end

