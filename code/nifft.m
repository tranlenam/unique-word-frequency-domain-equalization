function y=nifft(x)
y=ifft(x)*sqrt(length(x));