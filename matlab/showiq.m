function t=showiq(x,fs)

t=(0:length(x)-1)/fs;
plot(t,real(x)); hold on;
plot(t,imag(x));

hold off;