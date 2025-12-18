% Root Cosine Pulse Shaping
close all; clear; clc; fig_n = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare system parameters

% Symbol Generation:
% Consider using OMPSK here for better channel performance allegedly
numBits = 32;
numBytes = numBits / 8;
M = 4;              % M-PSK
numSymbols = numBits/log2(M);   % Number of symbols in transmission
Fs_in = 100; % symbol frequency (Hz)
symbol_duration = numSymbols/Fs_in; % input signal duration (s)

% IQ generation
carrier_frequency_multiplier = 1;                       % Carrier Sinusoid Cycles Per Symbol
carrier_frequency = carrier_frequency_multiplier*Fs_in; % Carrier Signal Frequency (Hz)
carrier_period = 1/carrier_frequency;                   % Carrier Signal Period (s)
sinusoid_sample_rate = 64;                             % Number of Samples Per Carrier Sinusoid
t_delta = Fs_in*carrier_frequency_multiplier*sinusoid_sample_rate; % Time Between Carrier Sinusoid Samples
carrier_to_input_clk_ratio = t_delta/Fs_in;          % Ratio of carrier clock to symbol Clock 
t = 0:1/t_delta:symbol_duration-1/t_delta;

% Pulse Shaping
L = 16;             % Polyphase Segments
span = 16;          % Filter span (pulse shape filter order is span*L)
rolloff = 1;      % Alpha (0 goes to sinc, 1 foes to more square shaped in time) "excess bandwidth"
F_int = 23.54;  % interpolation factor
Fs_out = Fs_in * F_int;  % output pulse shpaer frequency
NCO_bits = 0;
C = 2^NCO_bits;       % NCO modulus
phi_int = C/F_int;    % NCO Step

% Channel Parameters
linear_channel_attenuation = 1;   % greater number = greater attenuation
gauss_noise_level = .5;            % greater number = greater noise
salt_and_pepper_noise_level = .8;  % greater number = greater noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate MPSK symbols
bits = randi([0 1], numBits, 1);
% Plot the input bits and complex signals in time
symbol_axis = 0:1/Fs_in:symbol_duration*log2(M) - (1/Fs_in);
fig_n = stemplot(fig_n, symbol_axis, bits, 'Input BInary Signal', 'Time (s)', 'Magnitude');

% split into I and Q and assign constellation multiplier
[bitGroups, symbolMap, MpskSymbols] = genMpsk(M, bits);

% Plot the MPSK constellation 
fig_n = scatterplot(fig_n, real(symbolMap), imag(symbolMap), 'MPSK Constellation', 'In Phase Magntaude', 'Quadriture Magntaude', 'filled', 'MarkerFaceColor','b');

% Plot the input real and complex signals in time
symbol_axis = 0:1/Fs_in:symbol_duration - (1/Fs_in);
fig_n = stemplot(fig_n, symbol_axis, real(MpskSymbols), 'Input Mpsk I signal', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, imag(MpskSymbols), 'Input Mpsk Q signal', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate combined IQ signal
% If we want n sinusoids per symbol, we need to "sample and hold"
% the input symbol for n samples
real_up = repelem(real(MpskSymbols), carrier_to_input_clk_ratio);
complex_up = repelem(imag(MpskSymbols), carrier_to_input_clk_ratio);
symbol_axis_up = 0:1/t_delta:symbol_duration - (1/t_delta);

% Define I and Q
I = cos(2*pi*carrier_frequency*t);
Q = sin(2*pi*carrier_frequency*t);

% Multiply the input I and Q by the carrier I and Q
I_of_s = I .* real_up;
Q_of_s = Q .* complex_up;

% Plot the I and Q modulated signals
fig_n = stemplot(fig_n, t, I_of_s, 'I of s', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, t, Q_of_s, 'Q of S', 'Time (s)', 'Magnitude');

% Combine I and Q signals to form the modulated signal
s = I_of_s + Q_of_s;

% plot the summed modulated signal
fig_n = stemplot(fig_n, t, s, 's', 'Time (s)', 'Magnitude');

% Compute the Fourier Transform of the modulated signal
S = fft(s);
S = fftshift(S);
sy = (-length(S)/2:length(S)/2-1) * (Fs_in*carrier_to_input_clk_ratio/length(S));

% Plot the magnitude spectrum of the modulated signal
fig_n = stemplot(fig_n, sy, abs(S), 'Frequency Spectrum of IQ Mixer Output', 'Frequency (Hz)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate pulse shaping filter
h = genPulseFilter(rolloff, span, L, "sqrt");
[H, w] = freqz(h, 1, 4096);
H = [flip(H)' H']; w = [-flip(w)' w'];
H_dB = 20*log10(abs(H));

% Plot the impulse responses of the pulse shaping filter
fig_n = stemplot(fig_n, 1:length(h), h, 'RRC Impulse Response', 'Tap Index', 'Normalized Magnitude');

% Plot the frequency responses of the pulse shaping filter
fig_n = fig_n + 1;
figure(fig_n); % Full Pulse Shaping filter frequency response
plot(w, H_dB);
title('RRC Frequency Response');
set(gca, 'XTick', -pi:pi/10:pi);
set(gca, 'XTickLabel', {'-pi', '-0.9pi', '-0.8pi', '-0.7pi', '-0.6pi', ...
    '-0.5pi', '-0.4pi', '-0.3pi', '-0.2pi', '-0.1pi', '0', '0.1pi','0.2pi', ...
    '0.3pi', '0.4pi', '0.5pi', '0.6pi', '0.7pi', '0.8pi', '0.9pi', 'pi'});
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split pulse shaping filter into Polyphase Components and populate ROM
polyphase_coeffs = cell(1, L);
for k = 1:L
   polyphase_coeffs{k} = h(k:L:end);
end

row_size = length(h)/L;
column_size = L;

ROM = zeros(row_size, column_size);
for k = 1:column_size
    ROM(1:length(polyphase_coeffs{k}), k) = polyphase_coeffs{k}(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply pulse shaping filter to synthesised IQ output signal
shift_register = zeros(row_size,1);
output_signal = [];
accumulator = 0;

for k = 1:(length(s)+row_size-1) % ensure the whole signal runs through a zero initialized regester and runs all the way out
    if k <= length(s)
        shift_register = [real(s(k)); shift_register(1:end-1)];
    else
        shift_register = [0; shift_register(1:end-1)];
    end
    while accumulator < C
        ROM_addr_p = floor(L*accumulator/C);
        output_sum = sum(shift_register .* ROM(:,ROM_addr_p+1));
        output_signal = [output_signal, output_sum];
        accumulator = accumulator + phi_int;
    end
    accumulator = accumulator - C;
end
s_prime = output_signal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Outputs
t_delta_new = length(s_prime)/symbol_duration; % the signal has been interpolated by F_int
t_filtered = 0:1/t_delta_new:symbol_duration-1/t_delta_new;

% Plot the output of the pulse shaper filter
fig_n = stemplot(fig_n, t_filtered, s_prime, 'ROM-Shift Pulse Shaper Output', 'Time (s)', 'Magnitude');

% Compute the Fourier Transform of the output signal
S_prime = fft(s_prime);
S_prime = fftshift(S_prime);
f_S_prime = (-length(S_prime)/2:length(S_prime)/2-1) * (Fs_out*carrier_to_input_clk_ratio/length(S_prime));

% Plot the magnitude spectrum of the modulated signal
fig_n = stemplot(fig_n, f_S_prime, abs(S_prime), 'Frequency Spectrum of ROM-Shift Pulse Shaper Output', 'Frequency (Hz)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For DAC visualization purposes "convert the signal to analog"
% This assumes a DAC that doesnt do any interpolating
% Plot the output of the DAC
s_prime_t = s_prime;
fig_n = contplot(fig_n, t_filtered, s_prime_t, 'DAC Output', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% CHANNEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Gauss noise vector
gaussNoise = randn(length(s_prime),1)'/(1/gauss_noise_level);

% Generate Salt and Pepper noise vector
prob = 0.05; % probability of impulse
saltAndPepperNoise = zeros(1,length(s_prime));
impulses = rand(1,length(s_prime)) < prob;
saltAndPepperNoise(impulses) = (2*randi([0,1],1,sum(impulses))-1)'/(1/salt_and_pepper_noise_level);

total_noise = gaussNoise + saltAndPepperNoise;
fig_n = stemplot(fig_n, t_filtered, total_noise, 'Channel Noise', 'Time (s)', 'Magnitude');

s_prime_recieved_t = s_prime * (1/linear_channel_attenuation) + total_noise;

fig_n = contplot(fig_n, t_filtered, s_prime_recieved_t, 'Recieved Analog Channel Signal', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For ADC visualization purposes "convert the signal back to digital"
% Assume an ADC with a high sample rate
ADC_upsample_factor = 32;
% "sample a lot"
t_delta_ADC = t_delta_new*ADC_upsample_factor;
t_ADC = 0:1/t_delta_ADC:symbol_duration-1/t_delta_ADC;
% Check that the ADC has a LPF for interpolation. You may not get this for
% free
s_prime_recieved = interp(s_prime_recieved_t, ADC_upsample_factor);

% Plot the output of the ADC
fig_n = stemplot(fig_n, t_ADC, s_prime_recieved, 'ADC Output', 'Time (s)', 'Magnitude');

% Compute the Fourier Transform of the recieved signal
S_prime_recieved = fft(s_prime_recieved);
S_prime_recieved = fftshift(S_prime_recieved);
f_S_prime_recieved = (-length(S_prime_recieved)/2:length(S_prime_recieved)/2-1) * (Fs_out*carrier_to_input_clk_ratio*ADC_upsample_factor/length(S_prime_recieved));

% Plot the frequency spectrum of the recieved signal
fig_n = stemplot(fig_n, f_S_prime_recieved, abs(S_prime_recieved), 'Magnatude Spectrum of Recieved Signal', 'Frequency (Hz)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Defs
function [bitGroups, symbolMap, MpskSymbols] = genMpsk(M, bits)
    bitGroups = reshape(bits, log2(M), []).';
    symbolMap = genSymbolMap(M);
    if M==2
        MpskSymbols = (2*bitGroups - 1)';
    else 
        indices = bi2de(bitGroups, 'left-msb') + 1;
        grayIndices = bitxor(indices-1, floor((indices-1)/2)) + 1;
        MpskSymbols = symbolMap(grayIndices);
    end
end

function [symbolmap] = genSymbolMap(M)
    theta = pi/M;
    symbolmap = exp(1j*2*pi*(0:M-1)/M);
    if M~= 2
        symbolmap = symbolmap*exp(1j*theta);
    end
end

function [h] = genPulseFilter(rolloff, span, L, type)
    h = rcosdesign(rolloff, span, L, type);
    h = h(1:end-1); % trim one coeff so we have a rectangular ROM 
    h = h / max(h);
end

function fig_n = stemplot(fig_n, x, y, ttl, xlbl, ylbl, varargin)
    fig_n = fig_n + 1;
    figure(fig_n);
    stem(x,y,varargin{:});
    title(ttl);
    xlabel(xlbl); ylabel(ylbl);
    grid on;
end

function fig_n = contplot(fig_n, x, y, ttl, xlbl, ylbl, varargin)
    fig_n = fig_n + 1;
    figure(fig_n);
    plot(x,y,varargin{:});
    title(ttl);
    xlabel(xlbl); ylabel(ylbl);
    grid on;
end

function fig_n = scatterplot(fig_n, x, y, ttl, xlbl, ylbl, varargin)
    fig_n = fig_n + 1;
    figure(fig_n);
    scatter(x,y,varargin{:});
    title(ttl);
    xlabel(xlbl); ylabel(ylbl);
    grid on;
end