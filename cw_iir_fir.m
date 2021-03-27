%         Jacob Casey.  |  ECG Analysis.            % 
%          Start date 20/12/2020 20:23.             %

%Matlab cleaning
clc;
close all;
clear all;

% Read in files and set-up standard variables of the signal
filename = 'ecg.wav';
[y,Fs] = audioread(filename);
N = 42000;  
t = 0 : 1/Fs : (N - 1)*1/Fs; %create discrete time 

%% Apply IIR
    %set-up variables for IIR
    Ts = 1/Fs;
    fc = 0.1;       %cut-off frequency (Hz) 
    Fc = fc/Fs;     %Normalised cut-off freq
    order = 2;      % number of taps/order

    %Calc Co-eff
    [b,a] = butter(order, fc/(Fs/2),'high');


    %Apply filter to ecg
    y_IIR = filter(b,a,y);
    plot(y_IIR);
pause;
%% Apply FIR

%set-up variables

Ts = 1/Fs;              %sample period
fc = 10;                %cut off frequency in Hz 
Fc = fc/Fs;             %normalised cut off frequency
m = 20;                 %number of taps (N = 2m+1) 3.5/0.4
N = 2*m+1;              % total number of filter taps

for n = 1:m
    h(n) = 2*Fc*sin(n*2*pi*Fc)/(n*2*pi*Fc); %calculate truncated impulse response for LP filter (+ve n)
end

plot(h)
pause;
 
h = [fliplr(h) 2*Fc h];     %construct filter (add n = 0 coefficient for LP and -ve half)

figure;
w = hanning(N)';                    %generate N point hanning window
hw = h.*w;                          %apply window to filter coefficients
plot(hw)
title('Filter coefficients');
pause;

x = randn(1,1000000)*sqrt(512);     %generate white noise signal (normalise for default 512 point fft in pspectrum)
tic;
xf = conv(h,x);   %calculate filter output
toc;

tic;
hxf = conv(hw,x);
toc;
pause;
pspectrum(xf,Fs); hold; pspectrum(hxf,Fs); legend('Retangular Window','Hanning Window');

%Apply the FIR to the signal

tic;
ecg_IIR_FIR = conv(hw,y_IIR);
toc;

figure;
plot(ecg_IIR_FIR)
title('FIR filter');
pause;
