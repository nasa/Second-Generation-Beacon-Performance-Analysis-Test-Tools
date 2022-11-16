% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: chipscan 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: chipscan.m 11 2019-09-23 13:10:04Z reesebo $
% ///            

function [cro,v]=chipscan(id,fs,fchip,winsz,step,drw)

% rd = (sign(real(id)));
% id = (sign(imag(id)));
% indata=rd+1j*id;
indata=id;
% indata=upsample(indata(:).',2);
% fs=fs*2;

shift=round(fs/fchip/2);


bounds=[0.2*fchip, 10*fchip];
win=1:winsz;
for(ix=1:length(indata)/step-1)
    try
    cr(ix)=chiprate3(indata(win),fs,fchip,shift,bounds,0);
    win=win+step;
    t(ix) = win(fix(step/2)) * 1/fs;
    catch e
        break;
    end
end

% remove outliers
% for(ix=1:2)
%     u=mean(cr);
%     s=std(cr);
% 
%     idx=find(abs(cr-u) > 3*s);
% 
%     cr(idx)=[];
%     t(idx)=[];
% end

u=mean(cr); % update after last prune

% show results
if(drw)
    sfigure(drw);clf;hold on;
    plot(t,cr-fchip);
    title(['chip rate error over message length, u=' num2str(u)]);
    xlabel('time');
    ylabel('chiprate error (Hz)');
end

p1=polyfit(t,cr-fchip,1);
pv1=polyval(p1,t);
text(t(end),pv1(end),sprintf('y=%0.2fx+%0.2f', p1(1),p1(2)));
plot(t,pv1);
idx=find(t<0.167);
plot(t(idx),cr(idx)-fchip);
p2=polyfit(t(idx),cr(idx)-fchip,1);
pv2=polyval(p2,t(idx));
plot(t(idx),pv2);
text(t(idx(end)),pv2(end),sprintf('y=%0.2fx+%0.2f', p2(1),p2(2)));
cr1=mean(cr);
cr2=mean(cr(idx));

cro=[cr1 cr2];
v=[p1(1) p2(1)];
