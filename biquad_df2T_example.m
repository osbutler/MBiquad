clc; close all; clearvars;
% biquad_df2T_example - Example as for the biquad_df2T class. Stimulus is
% genereted and processed in blocks by the biquad_df2T class that uses 
% the specific biquad structure filter called "direct form 2 transposed".
% Before processing, the biquad filter type and parameters are defined.
% Output is then stored in a buffer. Last step is plotting modulus and 
% phase response of the filter and input and output.
% --------------------------
% Author:  Oscar Butler
% Project: MBiquad
% Date:    11.10.2023
% --------------------------
    
%% General settings
fs = 48e3;
BLOCKSIZE = 16; % Buffer size for ONE 32-bit floating point channel

%% Stimulus
N = 256;
Nfft = 2^14;
fSine = 1000;
t = 0:1/fs:(N-1)/fs;
f = 0:fs/Nfft:fs-1;
x = sin(2*pi*fSine*t);

%% Biquad parameters
lowshelf_param.type = 1;           % Lowshelf filter type
lowshelf_param.numStages = 1;       % Biquad order
lowshelf_param.gaindB = 0;         % Gain in dB
lowshelf_param.freqCut = 2500.0;    % Significant frequency 
lowshelf_param.Q = 0.707;          % Quality factor

%% Array init
inBuffer = zeros(1,BLOCKSIZE);

%% Init of LS filter instance
lowshelf = biquad_df2T(lowshelf_param, fs, BLOCKSIZE);
H = freqz(lowshelf.coeffs(1:3),[1 lowshelf.coeffs(4:5)],f,fs);

%% Block processing
for i=1:N/BLOCKSIZE-1
    inBuffer = x(i*BLOCKSIZE+1:(i+1)*BLOCKSIZE);
    lowshelf.mono_df2T(inBuffer);
    y(i*BLOCKSIZE+1:(i+1)*BLOCKSIZE) = lowshelf.outputBuffer;
end


%% Plot filter modulus and phase responses, stimulus and filtered stimulus
figure;
subplot(221)
semilogx(f,20*log10(abs(H)));
grid on;
xlim([10 20e3])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

subplot(222)
semilogx(f,180/pi*angle(H));
grid on;
xlim([10 20e3])
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')

subplot(2,2,[3 4])
plot(x)
hold on;
plot(y)
grid on;
xlabel('Samples')
ylabel('Amplitude')
xlim([0 N])

sgtitle('Filter frequency response and in/out signals in time domain')









