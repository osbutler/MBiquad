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
% Date:    11.4.2023
% --------------------------
    
%% General settings
fs = 48e3;
BLOCKSIZE = 128; % Buffer size for ONE 32-bit floating point channel

%% Stimulus
N = 2^12;
fSine = 50;
t = 0:1/fs:(N-1)/fs;
f = 0:fs/N:fs-1;
x = sin(2*pi*fSine*t);

%% Biquad parameters
type = 4;           % Lowshelf filter type
numStage = 1;       % Biquad order
gaindB = 3;         % Gain in dB
freqCut = 100.0;    % Significant frequency 
Q = 0.707;          % Quality factor

%% Array init
coeffLS = zeros(1,5);
stateLS = zeros(1,4*numStage);
inBuffer = zeros(1,BLOCKSIZE);
outBuffer = zeros(1,BLOCKSIZE);

%% Init of LS filter instance
LS = biquad_Stereo_df2T;
LS.init(numStage, coeffLS, stateLS, freqCut, Q, fs, gaindB, outBuffer, type);
LS.biquad_coeff_calculation;
H = freqz(LS.coeffs(1:3),[1 LS.coeffs(4:5)],f,fs);

%% Block processing
i = 1;
while i < N/BLOCKSIZE
    inBuffer = x(i*BLOCKSIZE+1:(i+1)*BLOCKSIZE);
    LS.mono_df2T(inBuffer, BLOCKSIZE);
    y(i*BLOCKSIZE+1:(i+1)*BLOCKSIZE) = LS.outputBuffer;
   i = i + 1; 
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
semilogx(f,angle(H));
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












