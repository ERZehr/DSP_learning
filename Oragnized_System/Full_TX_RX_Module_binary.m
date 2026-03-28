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
Fs_in = 1e7; % input symbol frequency (Hz)
BW_in = Fs_in / 2;
symbol_duration = numSymbols/Fs_in; % input signal duration (s)

% Pulse Shaping
L = 16;             % Polyphase Segments
span = 16;          % Filter span (pulse shape filter order is span*L)
rolloff = 1;      % Alpha (0 goes to sinc, 1 goes to more square shaped in time) "excess bandwidth"
F_int = 16;  %  pulse shpaer interpolation factor
Fs_out = Fs_in * F_int;  % output pulse shpaer sample frequency
BW_out = Fs_out / 2;
C = 4;
phi_int = C/F_int;    % NCO Step

% IQ generation
carrier_frequency = 2e8;

% Channel Parameters
linear_channel_attenuation = 1;   % greater number = greater attenuation
gauss_noise_level = 0.05;            % greater number = greater noise
salt_and_pepper_noise_level = 0;  % greater number = greater noise

% Farrow Parameters
Farrow_Length = 4;

% Sampler Parameters
rejection_ratio = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate bitstream and split into I and Q
bits = randi([0 1], numBits, 1);
% Plot the input bits and complex signals in time
symbol_axis = 0:1/Fs_in:numBits/Fs_in - (1/Fs_in);
fig_n = stemplot(fig_n, symbol_axis, bits, 'Input Binary Signal', 'Time (s)', 'Magnitude');

% split into I and Q and assign constellation multiplier
[bitGroups, symbolMap, MpskSymbols] = genMpsk(M, bits);

% Plot the input real and complex signals in time
symbol_axis = 0:1/Fs_in:symbol_duration - (1/Fs_in);
fig_n = stemplot(fig_n, symbol_axis, real(MpskSymbols), 'Input Mpsk I signal', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, imag(MpskSymbols), 'Input Mpsk Q signal', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and apply Barker interpolating pulse shaping filter
h = genPulseFilter(rolloff, span, L, "sqrt");
ROM = genPolyphaseROM(h, L, 'counter-clockwise');

% Apply pulse shaping filter to I
[I_symbols_shaped, Q_symbols_shaped] = applyBarkerInterpolator(ROM', real(MpskSymbols),imag(MpskSymbols),F_int,C);

% Plot the output of the pulse shaper filter
symbol_axis = 0:1/Fs_out:length(I_symbols_shaped)/Fs_out - (1/Fs_out);
fig_n = stemplot(fig_n, symbol_axis, I_symbols_shaped, 'I symbols shaped', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, Q_symbols_shaped, 'Q symbols shaped', 'Time (s)', 'Magnitude');

%DFTs
%fig_n = plotDFT(fig_n, I_symbols_shaped, Fs_out);
%fig_n = plotDFT(fig_n, Q_symbols_shaped, Fs_out);

% Eye Diagrams
% I
%figure(fig_n);
%eyediagram(I_symbols_shpaed, F_int);
%fig_n = fig_n + 1;
% Q
%figure(fig_n)
%eyediagram(Q_symbols_shpaed, F_int);
%fig_n = fig_n + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate combined IQ signal
[s_prime,~] = IQ_mixer(I_symbols_shaped,Q_symbols_shaped,carrier_frequency,Fs_out);

fig_n = stemplot(fig_n, symbol_axis, s_prime, 'Transmitted Signal', 'Time (s)', 'Magnitude');
% For DAC visualization purposes "convert the signal to analog"
fig_n = contplot(fig_n, symbol_axis, s_prime, 'DAC Output', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% CHANNEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Channel Noise
s_prime_recieved = apply_channel_noise(s_prime,gauss_noise_level,salt_and_pepper_noise_level,linear_channel_attenuation);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split I and Q  
[I_recovered,Q_recovered] = IQ_unmixer(s_prime_recieved,carrier_frequency,Fs_out);

fig_n = stemplot(fig_n, symbol_axis, I_recovered, 'Recieved I Channel Signal', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, Q_recovered, 'Recieved Q Channel Signal', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GOod up to here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply matched filter to I, Gardner TED, and Farrow filter
[I_symbols_recovered,Q_symbols_recovered, I_symbols_recovered_no_farrow, Q_symbols_recovered_no_farrow, ...
    TED_error_vector, TED_mu_values] = apply_timing_data_recovery(I_recovered,Q_recovered,ROM,F_int,Farrow_Length);

gardner_rate_axis = 0:1/Fs_out/F_int/2:length(I_symbols_recovered)/F_int/2/Fs_out - (1/Fs_out/F_int/2);
%fig_n = stemplot(fig_n, 1:length(TED_error_vector(1,:)), TED_error_vector(1,:), 'Gardner error vector I', 'Index', 'Magnitude');
%fig_n = contplot(fig_n, 1:length(TED_mu_values(1,:)), TED_mu_values(1,:), 'Gardner mu Values I', 'Index', 'Magnitude', '--');

fig_n = stemplot(fig_n, gardner_rate_axis, I_symbols_recovered, 'I Farrow Filter Output', 'Time (s)', 'Magnitude');
%fig_n = stemplot(fig_n, gardner_rate_axis, I_symbols_recovered_no_farrow, 'I Recovered NO Farrow Filter Output', 'Time (s)', 'Magnitude');
fig_n = fig_n + 1;
figure(fig_n);
plot(gardner_rate_axis, I_symbols_recovered, 'LineWidth', 1.5)
hold on
plot(gardner_rate_axis, I_symbols_recovered_no_farrow, 'LineWidth', 1.5)
hold off
title('I Values')
xlabel('Time')
ylabel('Amplitude')
legend('With Farrow', 'Without Farrow')
grid on

%fig_n = stemplot(fig_n, 1:length(TED_error_vector(2,:)), TED_error_vector(2,:), 'Gardner error vector Q', 'Index', 'Magnitude');
%fig_n = contplot(fig_n, 1:length(TED_mu_values(2,:)), TED_mu_values(2,:), 'Gardner mu Values Q', 'Index', 'Magnitude', '--');
fig_n = stemplot(fig_n, gardner_rate_axis, Q_symbols_recovered, 'Q Farrow Filter Output', 'Time (s)', 'Magnitude');
%fig_n = stemplot(fig_n, gardner_rate_axis, Q_symbols_recovered_no_farrow, 'Q Recovered NO Farrow Filter Output', 'Time (s)', 'Magnitude');
fig_n = fig_n + 1;
figure(fig_n);
plot(gardner_rate_axis, Q_symbols_recovered, 'LineWidth', 1.5)
hold on
plot(gardner_rate_axis, Q_symbols_recovered_no_farrow, 'LineWidth', 1.5)
hold off
title('Q Values')
xlabel('Time')
ylabel('Amplitude')
legend('With Farrow', 'Without Farrow')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample the recovered signal
h = gen_firpm_first_order_h('low', 0.49*pi, 0.51*pi, 1, 30);
ROM = genPolyphaseROM(h, 2, 'counter-clockwise');
[recieved_sampled_I, recieved_sampled_Q]= applyBarkerDecimator(ROM',I_symbols_recovered,Q_symbols_recovered,2);


recieved_sampled_symbol_axis_I = 0:1/Fs_in:length(recieved_sampled_I)/Fs_in - (1/Fs_in);
recieved_sampled_symbol_axis_Q = 0:1/Fs_in:length(recieved_sampled_Q)/Fs_in - (1/Fs_in);
fig_n = stemplot(fig_n, recieved_sampled_symbol_axis_I, recieved_sampled_I, 'Sampled I', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, recieved_sampled_symbol_axis_Q, recieved_sampled_Q, 'Sampled Q', 'Time (s)', 'Magnitude');

% implement a rejection ratio to reject transients 
recieved_sampled_I = recieved_sampled_I(abs(recieved_sampled_I) >= rejection_ratio);
recieved_sampled_Q = recieved_sampled_Q(abs(recieved_sampled_Q) >= rejection_ratio);

rejected_sampled_symbol_axis_I = 0:1/Fs_in:length(recieved_sampled_I)/Fs_in - (1/Fs_in);
rejected_sampled_symbol_axis_Q = 0:1/Fs_in:length(recieved_sampled_Q)/Fs_in - (1/Fs_in);
fig_n = stemplot(fig_n, rejected_sampled_symbol_axis_I, recieved_sampled_I, 'Sampled Accepted I', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, rejected_sampled_symbol_axis_Q, recieved_sampled_Q, 'Sampled Accepted Q', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% map the resulting symbols to data bits
recovered_symbols = recieved_sampled_I + 1j*recieved_sampled_Q;
recovered_bits = decMpsk(M, recovered_symbols, symbolMap);

recovered_symbols_axis = 0:1/Fs_in:length(recovered_symbols)/Fs_in - (1/Fs_in);
fig_n = stemplot(fig_n, 0:length(recovered_bits)-1, recovered_bits, 'Recovered bits', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BER = sum(bits ~= recovered_bits) / length(bits);
fprintf("Bit Error Rate: %f\n", BER);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%