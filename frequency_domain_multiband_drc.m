clc;clear all;close all;
Fs = 48000;
y = audioread('pureSpeech.wav');
y = stereo2mono(y);
% t = 1:10*48000;
% y = [sin(2*pi*500*t/Fs) 0.5*sin(2*pi*500*t/Fs) sin(2*pi*500*t/Fs)];
len=length(y);
frameSize=64;
hopSize=frameSize/4;
analysisWin=hamming(frameSize+1)'/1.08;
fftSize=frameSize;
analysisWin(end)=[];
ola=zeros(1,len);
ola1=zeros(1,len+fftSize);
numFrames=1+floor((len-frameSize)/hopSize);

% DRC Parameters. 
R = 2;
CT = -20;
att = hopSize/(Fs*0.2+1);
rel = hopSize/(Fs*0.5+1);
gain = 1;
log_gain = [];

% Subband Parameters.
numBinsPerBand = 8;
numBands = fftSize/2/numBinsPerBand;
bandBins = zeros(1,numBands);
BG = zeros(1,numBands);
VG = zeros(1,numBands);
binGains = zeros(1,fftSize/2+1);
for i=1:numFrames
    sig = y((i-1)*hopSize+1:(i-1)*hopSize+frameSize);
    sig_fft = rfft(sig,fftSize);
    magFFT = abs(sig_fft);
    for ind=1:numBands
       bandBins(ind) = sum(magFFT((ind-1)*numBinsPerBand+1:ind*numBinsPerBand));
    end
    bandBins(ind)=bandBins(ind)+magFFT(end);
    x_dB = 20*log10(sqrt(bandBins/fftSize));
    % WB-DRC gain application. 
    y_dB = CT + (x_dB-CT)/R;
    % ---------- Gain Calc and Smoothing ----------------
    BG = y_dB - x_dB;
    % Smoothing of gain based on attack and release times.
    for ind=1:numBands
        if(BG(ind)>VG(ind))
            VG(ind) = (1-att)*VG(ind) + att*BG(ind);
        else
            VG(ind) = (1-rel)*VG(ind) + rel*BG(ind);
        end
    end
    linGain = 10.^(VG/20);
    % Transform band gains to bin gains. 
    for ind=1:numBands
        binGains((ind-1)*numBinsPerBand+1:ind*numBinsPerBand) = linGain(ind);
    end
    binGains(end)=linGain(ind);
    sig_ifft = irfft(binGains.*sig_fft,fftSize);
    ola1((i-1)*hopSize+1:(i-1)*hopSize+fftSize) = ola1((i-1)*hopSize+1:(i-1)*hopSize+fftSize) + (hopSize/frameSize)*sig_ifft;
end
figure;plot(y); hold on; plot(ola1,'r');
% figure;plot(y);hold on; plot(ola1,'r');xlabel('Samples');ylabel('Original vs. OLA Synthesized Signal');
% figure;plot(y'-ola1(1:len));xlabel('Samples');ylabel('Error between Original and OLA Signal');

