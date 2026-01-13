% Root Cosine Pulse Shaping
close all; clear; clc; fig_n = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare system parameters

% Symbol Generation:
% Consider using OMPSK here for better channel performance allegedly
numBits = 64;
numBytes = numBits / 8;
M = 4;              % M-PSK
numSymbols = numBits/log2(M);   % Number of symbols in transmission
Fs_in = 1e8; % input symbol frequency (Hz)
symbol_duration = numSymbols/Fs_in; % input signal duration (s)

% Pulse Shaping
L = 16;             % Polyphase Segments
span = 16;          % Filter span (pulse shape filter order is span*L)
rolloff = 1;      % Alpha (0 goes to sinc, 1 goes to more square shaped in time) "excess bandwidth"
F_int = 16;  %  pulse shpaer interpolation factor
Fs_out = Fs_in * F_int;  % output pulse shpaer sample frequency
C = 1;
phi_int = C/F_int;    % NCO Step

% IQ generation
carrier_frequency = 1e8;

% Channel Parameters
linear_channel_attenuation = 1;   % greater number = greater attenuation
gauss_noise_level = .5;            % greater number = greater noise
salt_and_pepper_noise_level = .8;  % greater number = greater noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate MPSK symbols
bits = randi([0 1], numBits, 1);
% Plot the input bits and complex signals in time
symbol_axis = 0:1/Fs_in:numBits/Fs_in - (1/Fs_in);
fig_n = stemplot(fig_n, symbol_axis, bits, 'Input Binary Signal', 'Time (s)', 'Magnitude');

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
% Generate pulse shaping filter
h = genPulseFilter(rolloff, span, L, "normal");
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
% Apply pulse shaping filter to I and Q symbols
shift_register = zeros(row_size,1);
output_signal = [];
accumulator = 0;

I_symbols = real(MpskSymbols);

for k = 1:(length(I_symbols)+row_size-1) % ensure the whole signal runs through a zero initialized regester and runs all the way out
    if k <= length(I_symbols)
        shift_register = [real(I_symbols(k)); shift_register(1:end-1)];
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
I_symbols_shpaed = output_signal;


output_signal = [];
accumulator = 0;
Q_symbols = imag(MpskSymbols);

for k = 1:(length(Q_symbols)+row_size-1) % ensure the whole signal runs through a zero initialized regester and runs all the way out
    if k <= length(Q_symbols)
        shift_register = [real(Q_symbols(k)); shift_register(1:end-1)];
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
Q_symbols_shpaed = output_signal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the output of the pulse shaper filter
symbol_axis = 0:1/Fs_out:length(I_symbols_shpaed)/Fs_out - (1/Fs_out);
fig_n = stemplot(fig_n, symbol_axis, I_symbols_shpaed, 'I symbols shpaed', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, Q_symbols_shpaed, 'Q symbols shpaed', 'Time (s)', 'Magnitude');

% Compute the Fourier Transform of the I output signal
I_shpaed_dft = fft(I_symbols_shpaed);
I_shpaed_dft = fftshift(I_shpaed_dft);
I_shpaed_dft_axis = (-length(I_shpaed_dft)/2:length(I_shpaed_dft)/2-1) * (Fs_out/length(I_shpaed_dft));
% Compute the Fourier Transform of the Q output signal
Q_shpaed_dft = fft(Q_symbols_shpaed);
Q_shpaed_dft = fftshift(Q_shpaed_dft);
Q_shpaed_dft_axis = (-length(Q_shpaed_dft)/2:length(Q_shpaed_dft)/2-1) * (Fs_out/length(Q_shpaed_dft));

% Plot the magnitude spectrum of the modulated signals
fig_n = stemplot(fig_n, I_shpaed_dft_axis, abs(I_shpaed_dft), 'Frequency Spectrum of ROM-Shift Pulse Shaper Output I', 'Frequency (Hz)', 'Magnitude');
fig_n = stemplot(fig_n, Q_shpaed_dft_axis, abs(Q_shpaed_dft), 'Frequency Spectrum of ROM-Shift Pulse Shaper Output Q', 'Frequency (Hz)', 'Magnitude');

% Eye Diagrams
% I
figure(fig_n);
eyediagram(I_symbols_shpaed, F_int);
fig_n = fig_n + 1;
% Q
figure(fig_n)
eyediagram(Q_symbols_shpaed, F_int);
fig_n = fig_n + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate combined IQ signal
n = 0:length(I_symbols_shpaed)-1;

% Define I and Q
cos_lo = cos(2*pi*carrier_frequency/Fs_out*n);
sin_lo = sin(2*pi*carrier_frequency/Fs_out*n);
% Plot the LO signals
fig_n = stemplot(fig_n, symbol_axis, cos_lo, 'Cos LO', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, sin_lo, 'Sin LO', 'Time (s)', 'Magnitude');

% geneteate I and Q mixed Signals
I_mix = I_symbols_shpaed .* cos_lo - Q_symbols_shpaed .* sin_lo; % The actual real signal
Q_mix = I_symbols_shpaed .* sin_lo + Q_symbols_shpaed .* cos_lo;
% plot the I and Q mixed Signals
fig_n = stemplot(fig_n, symbol_axis, I_mix, 'I mix', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, Q_mix, 'Q mix', 'Time (s)', 'Magnitude');

s_prime = I_mix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For DAC visualization purposes "convert the signal to analog"
% This assumes a DAC that doesnt do any interpolating
% Plot the output of the DAC
s_prime_t = s_prime;
fig_n = contplot(fig_n, symbol_axis, s_prime_t, 'DAC Output', 'Time (s)', 'Magnitude');
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
fig_n = stemplot(fig_n, symbol_axis, total_noise, 'Channel Noise', 'Time (s)', 'Magnitude');

s_prime_recieved_t = s_prime * (1/linear_channel_attenuation) + total_noise;

fig_n = contplot(fig_n, symbol_axis, s_prime_recieved_t, 'Recieved Analog Channel Signal', 'Time (s)', 'Magnitude');
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