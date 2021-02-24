clc;clear all;close all;
Fs = 48000;
%y = audioread('pureSpeech.wav');
t = 1:10*48000;
y = [sin(2*pi*500*t/Fs) 0.5*sin(2*pi*500*t/Fs) sin(2*pi*500*t/Fs)];
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
VG = 0;
BG = 0;
log_gain = [];
for i=1:numFrames
    sig = y((i-1)*hopSize+1:(i-1)*hopSize+frameSize);
    sig_fft = fft(sig,fftSize);
    magFFT = abs(sig_fft);
    %rmsTime = 20*log10(sqrt(sum(sig.^2)/frameSize))
    rmsFreq = sqrt(sum(magFFT.^2))/fftSize;
    x_dB = 20*log10(sqrt(sum(magFFT.^2))/fftSize);
    % WB-DRC gain application. 
    y_dB = CT + (x_dB-CT)/R;
    % ---------- Gain Calc and Smoothing ----------------
    BG = y_dB - x_dB;
    % Smoothing of gain based on attack and release times.
    if(BG>VG)
        VG = (1-att)*VG + att*BG;
    else
        VG = (1-rel)*VG + rel*BG;
    end
    linGain = 10^(VG/20);
    log_gain = [log_gain linGain];
    sig_ifft = ifft(linGain*sig_fft,fftSize);
    ola1((i-1)*hopSize+1:(i-1)*hopSize+fftSize) = ola1((i-1)*hopSize+1:(i-1)*hopSize+fftSize) + (hopSize/frameSize)*sig_ifft;
end
figure;plot(y); hold on; plot(ola1,'r');
% figure;plot(y);hold on; plot(ola1,'r');xlabel('Samples');ylabel('Original vs. OLA Synthesized Signal');
% figure;plot(y'-ola1(1:len));xlabel('Samples');ylabel('Error between Original and OLA Signal');

