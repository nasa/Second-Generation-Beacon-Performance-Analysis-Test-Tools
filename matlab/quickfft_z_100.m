function [ X1, faxis ] = quickfft_z(varargin)
% Function:	quickfft()
%
% Purpose:	Computes and displays the fft of the array x
%
% Syntax:	X = quickfft(x,srate);
%
% Arguments:	srate:		Sample rate
%		x:		Input data array
%
% Output:	X:		FFT(x)
%
% Notes:	1)  routine computes fft using power of two closest
%		    to length of input array up to 8192.
x=varargin{1};
srate=varargin{2};
if(nargin==2)
    plt='b';
else
    plt = varargin{3};
end 

n = length(x);
N = 2.^(fix(log10(n)/log10(2)));
N = fix(min(N,131072)/100);
if(nargin > 3)
        N=varargin{4};  
end
if(nargin > 4)
    cf=varargin{5};
end

if(1)
W = blackman(N);
X = fft(x(1:N).*W');
else
    W = window(@blackmanharris,N);
    X = fft(x(1:N).*W.');
end
X1 = 20*log10(abs(fftshift(X))+eps)-10*log10(srate/N);
%X1=fliplr(X1);
%X1=X1-max(X1)-10;
faxis = [-N/2:N/2-1]*srate/N;
if(exist('cf') == 1)
    faxis = (faxis+cf)/1e3;
else
    faxis = (faxis)/1e3;
end
stdX = std(X1);
meanX = mean(X1);
V = [ faxis(1) faxis(N) (meanX-2*stdX) max(X1) ];
axis(V);
plot(faxis,medfilt1(X1,2),'linewidth',1);
grid;
xlabel('Frequency KHz');
ylabel('Power - dB');
axis;
drawnow;
