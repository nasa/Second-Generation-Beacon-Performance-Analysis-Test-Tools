% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: pk2pk 
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 11 $ $Date: 2019-09-23 09:10:04 -0400 (Mon, 23 Sep 2019) $
% ///				$Id: pk2pk.m 11 2019-09-23 13:10:04Z reesebo $
% ///            
function pk2pko=pk2pk(xoc,fs,fchip,fig)

    it=round(fs/fchip/2);
    for(os=1:it)
        [pr,pi]=getpk(xoc,fs,fchip,os,0);
        idx=find(pr>0);
        pt=pr(idx);
        km = kmeans(pt(:),2);
        idx1=find(km==1);  % find cluster of greater value
        idx2=find(km==2);
        dx(os) = abs(mean(pr(idx(idx1))) - mean(pr(idx(idx2))));
    end
    [~,os]=min(dx);
    
    [pr,pi]=getpk(xoc,fs,fchip,os,fig);
    
    mx=mean(pr(pr>0));
    mn=mean(pr(pr<0));
    dr = abs(mx-mn);
  
    mx=mean(pi(pi>0));
    mn=mean(pi(pi<0));
    di = abs(mx-mn);

    pk2pko = ((dr/di-1)*100);
    
end

function [pr,pi]=getpk(xoc,fs,fchip,os,fig)


re=real(xoc(1:fs));
im=imag(xoc(1:fs));

it = fs/fchip;

pr=getpint(re,it,os,fs,fig);
pi=getpint(im,it,os+it/2,fs,fig);
end

function pr=getpint(re,it,os,fs,fig)
win=(1:it)+os;

pr=zeros(length(fs/4),1);
try
for(ix=1:fs/it)
    pr(ix)=sum(re(win))/4;
    win=win+it;
end
catch
end

t1=(0:fs-1)/fs;
t2=((0:4:fs-1)+os+it/2)/fs;
if(fig)
    %sfigure(fig);clf;hold on;
    clf;hold on;
    plot(t1,re);
    plot(t2(1:length(pr)),pr,'o');
    axis([0.036246123602067   0.037447855039413  -0.954393075519784   0.874565264162429]);
end
disp('');
end