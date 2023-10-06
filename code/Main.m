clear all
clc
M=4;        % QPSK
k=log2(M);  % # of bits/symbol

% iid Rayleigh channel
channellength=15;
powerprofile=(1/channellength)*ones(1,channellength); % equal power/tap
powerprofile=sqrt((powerprofile.')/2); % normalization


blocklength=64;
uwlength=16;
datalength=blocklength-uwlength;


nIterations=1e4;

EbNo=[0:2:14];
EsNo=EbNo+10*log10(k); % dB scale
EsNo_eff=EsNo-10*log10(1+uwlength/(blocklength-uwlength));
Es=1;
noisepower=sqrt(Es./(10.^(EsNo_eff/10))/2); % compute the noise power

uniqueword=reshape((randn(uwlength*k,1)>0),[],k); % Gaussian distribution
uniqueword=num2str(uniqueword);
uniqueword=bin2dec(uniqueword);
uniqueword=pskmod(uniqueword,M); % modulate data

totallength=blocklength+uwlength+channellength;


BER=zeros(length(EbNo),1);
for iEbNo=1:length(EbNo)
    nErrors = 0;
    for iIteration=1:nIterations
        % 0/1 data
        bindata=(randn(datalength*k,1)>0);
        bindata=double(bindata);
        decdata=reshape(bindata,[],k);
        decdata=(bit2int(decdata',2))';
        moddata=pskmod(decdata,M,pi/4);
        % Form the data block to be transmitted
        datablock=[uniqueword;moddata;uniqueword];
        
        channel=powerprofile.*(randn(channellength,1)+1i*randn(channellength,1)); 
        
        % compute the channel output and add Gausian noise
        rxdata=conv(datablock,channel)+noisepower(iEbNo)*(randn(totallength-1,1)...
           +1i*randn(totallength-1,1));
        rxdata=rxdata(uwlength+1:uwlength+blocklength);
        
        
        % Receiver
        z =optimalHDFEnoisepredictor(rxdata,uniqueword,M,channel,2*noisepower(iEbNo)^2);
        % Demodulate data
        demoddata=pskdemod(z,M,pi/4);
        demoddata=demoddata(1:length(decdata));
        
        bindetectdata=(int2bit(demoddata',2))';
        bindetectdata=reshape(bindetectdata,[],1);
        nErrors=nErrors+sum(bindata~=bindetectdata);
    end
    BER(iEbNo)=nErrors/(length(bindata)*nIterations);
end
semilogy(EbNo,BER)
xlabel('EbNo')
ylabel('BER')
saveas(gcf,'../results/BER.png')