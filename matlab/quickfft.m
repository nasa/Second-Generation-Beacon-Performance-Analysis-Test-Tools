function [ X1, faxis, X ] = quickfft(varargin)
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
    plt=[];
    col=[];
elseif(nargin==4)
    plt = varargin{3};
    col=varargin{4};
else
    plt = varargin{3};
    col=[];
end

sz=size(x);
if(sz(1) > sz(2))
    x=x.';
end
n = length(x);
N = 2.^(fix(log10(n)/log10(2)));

%N = min(N,131072);
W = blackman(N);
X = fft(x(1:N).*W'/sum(W));
X1 = 20*log10(abs(fftshift(X))+eps);
faxis = [-N/2:N/2-1]*srate/N/1e3;
stdX = std(X1);
meanX = mean(X1);
V = [ faxis(1) faxis(N) (meanX-2*stdX) max(X1) ];
try
%axis(V+eps);
catch
end
% if(col)
%     %plot(faxis,X1,plt,'color',col);
%     plot(faxis,X1);
% else
%     plot(faxis,X1,plt);
% end
%plot(faxis,medfilt1(X1,20)); %,'linewidth',3);
if(~isempty(plt))
    plot(faxis,X1,plt); %,'linewidth',3);
else
    plot(faxis,X1); %,'linewidth',3);
end
grid;
xlabel('Frequency (kHz)');
ylabel('Power dB');
axis;
drawnow;
