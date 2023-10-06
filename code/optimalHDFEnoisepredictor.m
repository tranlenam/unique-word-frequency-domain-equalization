function equalizedsignal = optimalHDFEnoisepredictor(rxdata,uniqueword,...
    constellationsize,channel,noisepower)
uwlength=length(uniqueword);
blocklength=length(rxdata);

% Fourier Transform;
Iprime=[eye(blocklength-uwlength) zeros(blocklength-uwlength,uwlength);...
       zeros(uwlength,blocklength-uwlength) zeros(uwlength)];

Omega=nfft((nfft(Iprime))')'; % The matrix Omega defined below (7)

U=[zeros(blocklength-uwlength,1);uniqueword];
U=nfft(U); % normalized FFT


channellength=length(channel);
fftofchannel=fft(channel,blocklength);
D=diag(fftofchannel); 

G=Omega*D'*inv(D*Omega*D'+noisepower*eye(blocklength)); % matrix G defined below (7)

Psi=diag(Omega-G*D*Omega); % The matrix Psi defined below (8)

% Nb =2;   % # of feedback taps.
% Nb =5;   % # of feedback taps.
Nb = channellength-1;   % # of feedback taps.
v=blocklength*ifft(Psi,blocklength); 
lambda=v(2:1+Nb);   % vector lambda in (10)
Lambda=toeplitz(v(1:Nb)');
b=Lambda\lambda;


equalizedsignal=zeros(blocklength,1)+1i*zeros(blocklength,1);
% Feedbach Register 
reg=zeros(Nb,1);

% Initialize the feedback register with the Unique word

fftofrxdata = 1/sqrt(blocklength)*fft(rxdata);
fftofrxdata = fftofrxdata-D*U; % cancel the interference caused by UW

% Transform to time domain signal
r=sqrt(blocklength)*ifft(G*fftofrxdata);
reg(1:Nb)=r(end:-1:(end-Nb+1));
%*************************************************************************%

for k=1:length(rxdata)
    equalizedsignal(k)   =r(k)-sum(b.*reg);
    detectedsymbol   =pskdemod(equalizedsignal(k),constellationsize,pi/4);
    detectedsymbol   =pskmod(detectedsymbol,constellationsize,pi/4);
    % shift the register;
    reg                     =circshift(reg,1);
    reg(1)                  =r(k)-detectedsymbol;            
end
equalizedsignal = equalizedsignal(1:end-uwlength);