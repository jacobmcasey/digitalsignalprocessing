%         Jacob Casey.  |  ECG Analysis.            % 
%          Start date 20/12/2020 20:23.             %

%Matlab cleaning
clc;
close all;
clear all;

% Read in files and set-up standard variables of the signal
filename = 'ecg.wav';
[y,Fs] = audioread(filename);
N = Fs; %set N to Fs to achieve 1hz bin resolution

pspectrum(y);

%Apply FFT with N bins
fft_y = fft(y,N);

%Remove the repeated half
fft_y = fft_y(1:length(fft_y)/2+1);

%Calculate power
power = 1/(Fs*N) * abs(fft_y).^2;

%Set k so we can plot FFT from zero
k=0:1:length(fft_y)-1;

%Plot with logarithmic scale
plot(10*log10(power));

nn = hanning(N/2+1);       
mm = hamming(N/2+1); %generate N point hanning window

nnw = power.*nn;
mmw = power.*mm;

figure;
plot(k,10*log10(nnw)); hold;
plot(k,10*log10(mmw)); legend('hanning window','hamming window')
title('Implementing FFT using own code')
xlabel('Frequency (Hz)')
ylabel('Power Spectrum(Db)')  ;

