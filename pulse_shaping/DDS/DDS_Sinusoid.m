%% Sampling Practice
clear; clc; close all; fig_n = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clock Speed
clk = 400e6;
% pre generated
lut_amplitude_resolution        = 16; % bit resolution per sample in the ROM
lut_period_resolution           = 14; % bits in period index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create sinusoidal ROM (This would be an array LUT at compile time)
lut_period_size = 2^lut_period_resolution; % number of samples per period in the ROM
ROM_bytes       = lut_amplitude_resolution*lut_period_size/8; % number of bytes the ROM takes up

% Generate the ROM
n        = 0:lut_period_size-1;
sine_rom = sin(2*pi*n/lut_period_size);

% Plot the sinusoidal LUT
fig_n = stemplot(fig_n, n, sine_rom, 'Sine ROM', 'Sample Index', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
frequency_multiplier = 3.14; % sent in DDS request packet
frequency_multiplier_slope = 0.000001; % sent in DDS request packet

initial_phase = .5; % sent in DDS request packet

amplitude_multiplier = 1; % sent in DDS request packet
amplitude_multiplier_slope = 0.0001; % sent in DDS request packet

stride = 1; % Define the duration for the generation process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DDS Parameters
accumulator = initial_phase;
phi_int = frequency_multiplier/lut_period_size;

output_frequency = zeros(1,floor(lut_period_size*stride));
modulated_signal = zeros(1,floor(lut_period_size*stride));
%%
% Modulate the input signal with the NCO signal
for k = 1:floor(lut_period_size*stride)
    if accumulator >= 1
        accumulator = accumulator - 1;
    else
        ROM_index = mod(floor(accumulator*lut_period_size), lut_period_size) + 1;
        modulated_signal(k) = amplitude_multiplier * sine_rom(ROM_index);
        output_frequency(k) = frequency_multiplier*clk/lut_period_size;

        amplitude_multiplier = amplitude_multiplier + amplitude_multiplier_slope;
        phi_int = phi_int + frequency_multiplier_slope;
    end
    accumulator = accumulator + phi_int;
end

% Plot the modulated signal and its frequency
t = 0:1/clk:length(modulated_signal)/clk-1/clk;
fig_n = stemplot(fig_n, t, modulated_signal, 'DDS Signal', 'Time (s)', 'Magnitude');
fig_n = contplot(fig_n, t, output_frequency, 'DDS Frequency', 'Time (s)', 'Frequency');

% Compute the Fourier Transform of the modulated signal
modulated_signal_dft = fft(double(modulated_signal));
modulated_signal_dft = fftshift(modulated_signal_dft) / length(modulated_signal);
freq_w = (-length(modulated_signal_dft)/2:length(modulated_signal_dft)/2-1) .* frequency_multiplier*clk/length(modulated_signal_dft);

% Plot the magnitude spectrum of the modulated signal
fig_n = contplot(fig_n, freq_w, abs(modulated_signal_dft), 'Magnitude Spectrum of DDS Signal', 'Frequency (Hz)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Defs
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