function y=nfft(x)
y=fft(x)./sqrt(length(x));