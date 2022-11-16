% ///	Copyright(c) 2017 United States Government as represented by the 
% ///	Administrator for The National Aeronautics and Space Administration.  
% ///	All Rights Reserved. 
% ///	
% ///		Government Agency: NASA 
% ///		Government Agency Original Software Designation: GSC-18375-1
% ///		Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
% ///		User Registration Requested.  Please Visit https://software.nasa.gov/
% ///     
% ///     Module: correcttilt.m
% ///     
% ///     Author:   Reese Bovard
% ///             Concentric Real Time, LLC
% ///   
% ///     [version]:	$Revision: 15 $ $Date: 2022-09-29 11:45:13 -0400 (Thu, 29 Sep 2022) $
% ///				$Id: correcttilt.m 15 2022-09-29 15:45:13Z reesebo $
% ///            

function rxx=correcttilt(rx,fig)
len=length(rx);
rxx=zeros(1,len);
center=fix(0.5*len); %fix(0.1*len);
win=1:fix(center/4);
win=win+fix(center/2);
% win=10000:20000;
re=real(rx(win));
im=imag(rx(win));
cx=mean(re(re>0&im>0));
cy=mean(im(re>0&im>0));
tilt=pi/4-atan2(cx,cy);
rxx = rx .* exp(-1j*tilt);

fprintf('tilt=%d\n', tilt)

if(fig)
    sfigure(fig);clf;hold on;
    
    plot(rx,'.');
    plot(rxx,'.');
    plot(cx,cy,'x');
    Constellation=.7071*[1+1i,1-1i,-1+1i,-1-1i];

    plot(Constellation,'gd');
    title('symbol map center to end');
    axis([-1 1 -1 1]);
    grid;
end