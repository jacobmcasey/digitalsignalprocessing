% Jacob Casey. 2020. 
% Generate 2 sinusoids then FFT to find frequencies with a FFT resolution
% of 1 Hz 
%
% Hard coded funtion to get DFT matrix, then timed to compare speed difference
% between FFT and DFT.
%
%

j = sqrt(-1);

fs = 32;            %sampling frequency
Ts = 1/fs;          %sample period
f1 = 8;            %frequency of 1st sinusoid
a1 = 1;           %amplitude of 1st sinusoid
f2 = 10;            %frequency of 1st sinusoid
a2 = 1;           %amplitude of 1st sinusoid
N = 32;             %no of samples in transform

t = 0:Ts:(N-1)*Ts; 
%Rectangle FFT function 

s1 = a1*cos(2*pi*f1*t) + a2*cos(2*pi*f2*t);


figure; plot(t,s1);    %plot
xlabel('t'); ylabel('s1'); title('Discrete time sinusoid (f1 + f2)');
grid on

tic
s1_fft = fft(s1);
toc
k = 0:1:N-1;

fft_real = real(s1_fft);

figure;
plot(k,abs(s1_fft));     %plot magnitude of FFT coefficients

sum = 0.0;
s1_dft = zeros(1,N);

tic
for K = 0:N-1
    sum = 0;
    
    for n = 0:N-1
        W = exp((-j*2*pi*K*n)/N);
        sum = sum + s1(n+1) * W;
    end
    
    s1_dft(K+1) = sum;

end
toc

figure;
plot(k,abs(s1_dft)); 


