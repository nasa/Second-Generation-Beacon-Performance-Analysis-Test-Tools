% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: evalchiperror 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: evalchiperror.m 15 2022-09-29 15:45:13Z reesebo $
% ///            
function [cro,v]=evalchiperror(err,fs,fchip,win,fig)

err=err(:).';

fac=fs/fchip;
err=unwrap(err*2*pi)/pi/2/fac; % factor of 1/4 normalizes to chips given fs is 4x chips.
samp=fix(win*fs);
% win=1:samp;
% erro=[];
% try
% for(ix=1:length(err)/samp*2)
%     erro(ix) = mean(err(win));
%     win=win+samp/2;
% end
% catch
% end
% err=erro;
% fs=fs/samp*2;
err=err-mean(err);
tx=1/fs*(0:length(err)-1);tx=tx(:).';
if(fig)
    sfigure(fig);clf;hold on;
    plot(tx,err);
end
% %plot(tx(2:end),diff(error)/4);
% p=polyfit(tx(:),err(:),2);
% plot(tx,polyval(p,tx));

[p1,s,mu]=polyfit(tx,err,1);
cro(1) = p1(1)/mu(2);
pv1=polyval(p1,tx,s,mu);
idx=find(tx<1 & tx>0);

[p1_2,s,mu]=polyfit(tx(idx),err(idx),2);
v(1)=p1_2(1)/mu(2)^2 *2;% + p1_2(2)/mu(2); % *(tx(idx(end))-mu(1));
pv1=polyval(p1_2,tx(idx),s,mu);
%plot(tx(idx),pv1);

if(fig)
    text(tx(end),pv1(end),sprintf('e=%0.2f,v=%0.2f', cro(1),v(1)));
    plot(tx(idx),pv1);
end
idx=find(tx<0.167 & tx>0.07);
plot(tx(idx),err(idx));
[p2,s,mu]=polyfit(tx(idx),err(idx),1);
cro(2) = p2(1)/mu(2);
pv2=polyval(p2,tx(idx),s,mu);

[p2_2,s,mu]=polyfit(tx(idx),err(idx),2);
v(2)=p2_2(1)/mu(2)^2 *2; % * (tx(idx(end))-mu(1));

plot(tx(idx),pv2);
text(tx(idx(end)),pv2(end),sprintf('e=%0.2f,v=%0.2f', cro(2),v(2)));

cro=[p1(1) p2(1)];
v=[p1_2(1)*2 p2_2(1)*2];
if(fig)
    title(sprintf('Chip Rate Error: Timing Error Detector'));
    xlabel('Time (s)');
    ylabel('Error (chips)');grid on;
    legend('error','rate');
end
