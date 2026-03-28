% 1 comms channel, 7 RF channels
close all; clear; clc; fig_n = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare system parameters

% Signal Parameters
Fs_in = 5000; % input symbol frequency (Hz)
t = 0:1/Fs_in:1-1/Fs_in;


% Comms Symbol Generation:
% Consider using OMPSK here for better channel performance allegedly
numBits = 4096;
numBytes = numBits / 8;
M = 4;              % M-PSK
numSymbols = numBits/log2(M);   % Number of symbols in transmission
BW_in = Fs_in / 2;
symbol_duration = numSymbols/Fs_in; % input signal duration (s)
rejection_ratio = 0.1;


% Pulse Shaping/Barker Interpolator Parameters
L = 8;             % Polyphase Segments
span = 8;          % Filter span (pulse shape filter order is span*L)
rolloff = 0.8;     % Alpha (0 goes to sinc, 1 goes to more square shaped in time) "excess bandwidth"
F_int = 16;        %  pulse shpaer interpolation factor
Fs_out = Fs_in * F_int;  % output pulse shpaer sample frequency
BW_out = Fs_out / 2;
C = 4;
phi_int = C/F_int;    % NCO Step


% Channelizer Parameters
Fs_in = 100000; % 0-24kHz range
K = 8; % number of channels
D = 8; % number of decimations
% Prototype Filter
wp = pi/K-pi/100;
ws = pi/K+pi/100;
rp = 0.1;
rs = 100;


% Channel Parameters
linear_channel_attenuation = 1;   % greater number = greater attenuation
gauss_noise_level = 0;            % greater number = greater noise
salt_and_pepper_noise_level = 0;  % greater number = greater noise


% Farrow Parameters
Farrow_Length = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare Input Signals
x0 = randi([0 1], numBits, 1);
x1 = cos(2*pi*100.*t) + cos(2*pi*200.*t);
x2 = cos(2*pi*100.*t) + sin(2*pi*200.*t);
x3 = sin(2*pi*100.*t) + cos(2*pi*200.*t);
x4 = sin(2*pi*100.*t) + sin(2*pi*200.*t);
x5 = cos(2*pi*100.*t) + 2*cos(2*pi*200.*t);
x6 = cos(2*pi*100.*t) + 2*sin(2*pi*200.*t);
x7 = sin(2*pi*100.*t) + 2*cos(2*pi*200.*t);

channels = {x0, x1, x2, x3, x4, x5, x6, x7};

for channel = 1:8
    fig_n = stemplot(fig_n, 1:length(channels{channel}), channels{channel}, ['Input Signal x' num2str(channel-1)], 'Index', 'Magnitude');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put Digital Data throught the Constellation Mapper
[bitGroups, symbolMap, x0_symbols] = genMpsk(M, x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and apply Barker interpolating pulse shaping filter
h = genPulseFilter(rolloff, span, L, "sqrt");
Interpolation_ROM = genPolyphaseROM(h, L, 'counter-clockwise');

% Apply pulse shaping filter to comms data
[comms_I_symbols_shaped, comms_Q_symbols_shaped] = applyBarkerInterpolator(Interpolation_ROM', real(x0_symbols),imag(x0_symbols),F_int,C);
% Apply pulse shaping filter to RF data
[RF1_symbols_shaped, ~] = applyBarkerInterpolator(Interpolation_ROM', x1,zeros(1,length(x1)),F_int,C);
[RF2_symbols_shaped, ~] = applyBarkerInterpolator(Interpolation_ROM', x2,zeros(1,length(x2)),F_int,C);
[RF3_symbols_shaped, ~] = applyBarkerInterpolator(Interpolation_ROM', x3,zeros(1,length(x3)),F_int,C);
[RF4_symbols_shaped, ~] = applyBarkerInterpolator(Interpolation_ROM', x4,zeros(1,length(x4)),F_int,C);
[RF5_symbols_shaped, ~] = applyBarkerInterpolator(Interpolation_ROM', x5,zeros(1,length(x5)),F_int,C);
[RF6_symbols_shaped, ~] = applyBarkerInterpolator(Interpolation_ROM', x6,zeros(1,length(x6)),F_int,C);
[RF7_symbols_shaped, ~] = applyBarkerInterpolator(Interpolation_ROM', x7,zeros(1,length(x7)),F_int,C);

% Plot outputs
shaped_comms_channels = {comms_I_symbols_shaped, comms_Q_symbols_shaped};
shaped_RF_channels = {RF1_symbols_shaped, RF2_symbols_shaped, RF3_symbols_shaped, RF4_symbols_shaped, RF5_symbols_shaped, RF6_symbols_shaped, RF7_symbols_shaped};

fig_n = stemplot(fig_n, 1:length(shaped_comms_channels{1}), shaped_comms_channels{1}, 'Input Comms Signal I', 'Index', 'Magnitude');
fig_n = stemplot(fig_n, 1:length(shaped_comms_channels{2}), shaped_comms_channels{2}, 'Input Comms Signal Q', 'Index', 'Magnitude');
for channel = 1:7
    fig_n = stemplot(fig_n, 1:length(shaped_RF_channels{channel}), shaped_RF_channels{channel}, ['Input RF Signal x' num2str(channel)], 'Index', 'Magnitude');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate combined Comms IQ signal
[mixed_comms_signal,~] = IQ_mixer(comms_I_symbols_shaped,comms_Q_symbols_shaped,Fs_in,Fs_out);
fig_n = stemplot(fig_n, 1:length(mixed_comms_signal), mixed_comms_signal, 'Mixed Comms SIgnal', 'Index', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Synthesis Channelizer
% Generate prototype filter
h0 = gen_firpm_first_order_h('low', wp, ws, rp, rs);
h0 = [h0 zeros(1, mod(-length(h0), K))];
ROM = genPolyphaseROM(h0, K, 'counter-clockwise');

% Pad signal vectors out to the same length 
signals_to_inv_channelizer = {mixed_comms_signal; RF1_symbols_shaped; RF2_symbols_shaped; RF3_symbols_shaped; RF4_symbols_shaped; RF5_symbols_shaped; RF6_symbols_shaped; RF7_symbols_shaped};
max_length = max(cellfun(@length, signals_to_inv_channelizer));
signals_padded = cellfun(@(v) [v(:).', zeros(1, max_length - length(v))], signals_to_inv_channelizer, 'UniformOutput', false);
signals_to_inv_channelizer = vertcat(signals_padded{:});


[merged_signals_to_transmit, ~] = applyPolyphaseIFFTInverseChannelizer(ROM, signals_to_inv_channelizer, zeros(size(signals_to_inv_channelizer)), K);
fig_n = stemplot(fig_n, 1:length(merged_signals_to_transmit), merged_signals_to_transmit, 'All Merged Signals Discrete', 'Index', 'Magnitude');
fig_n = contplot(fig_n, 1:length(merged_signals_to_transmit), merged_signals_to_transmit, 'All Merged Signals Analog', 'Index', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%% CHANNEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Channel Noise
recieved_signal = apply_channel_noise(merged_signals_to_transmit, gauss_noise_level, salt_and_pepper_noise_level, linear_channel_attenuation);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Analysis Channelizer
ROM = genPolyphaseROM(h0, K, 'clockwise');
[analysis_channelizer_output, ~] = applyPolyphaseFFTChannelizer(ROM, recieved_signal, zeros(1,length(recieved_signal)), K);
analysis_channelizer_output = K * analysis_channelizer_output;

fig_n = stemplot(fig_n, 1:length(analysis_channelizer_output(1,:)), analysis_channelizer_output(1,:), 'Analysis Channelizer Comms Channel Mixed', 'Index', 'Magnitude');
for channel = 2:8
    fig_n = stemplot(fig_n, 1:length(analysis_channelizer_output(channel,:)), analysis_channelizer_output(channel,:), ['Analysis Channelizer RF Signal x' num2str(channel - 1)], 'Index', 'Magnitude');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split I and Q  
[Recieved_Comms_Signal_I,Recieved_Comms_Signal_Q] = IQ_unmixer(analysis_channelizer_output(1,:),Fs_in,Fs_out);

fig_n = stemplot(fig_n, 1:length(Recieved_Comms_Signal_I), Recieved_Comms_Signal_I, 'Recieved I Comms Signal', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, 1:length(Recieved_Comms_Signal_Q), Recieved_Comms_Signal_Q, 'Recieved Q Comms Signal', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conduct Timing and Error Detection
ROM_matched = Interpolation_ROM;
[Comms_I_symbols_recovered,Comms_Q_symbols_recovered, Comms_I_symbols_recovered_no_farrow, Comms_Q_symbols_recovered_no_farrow, Comms_TED_error_vector, Comms_TED_mu_values] = apply_timing_data_recovery(Recieved_Comms_Signal_I,Recieved_Comms_Signal_Q,ROM_matched,F_int,Farrow_Length);

RF_symbols_recovered_vec = [];
RF_symbols_recovered_no_farrow_vec = [];
RF_TED_error_vector_vec = [];
RF_TED_mu_values_vec = [];

for channel = 2:8
    [RF_symbols_recovered,~, RF_symbols_recovered_no_farrow, ~, RF_TED_error_vector, RF_TED_mu_values] = apply_timing_data_recovery(analysis_channelizer_output(channel,:),zeros(1,length(analysis_channelizer_output(channel,:))),ROM_matched,F_int,Farrow_Length);
    RF_symbols_recovered_vec = [RF_symbols_recovered_vec; RF_symbols_recovered];
    RF_symbols_recovered_no_farrow_vec = [RF_symbols_recovered_no_farrow_vec; RF_symbols_recovered_no_farrow];
    RF_TED_error_vector_vec = [RF_TED_error_vector_vec; RF_TED_error_vector];
    RF_TED_mu_values_vec = [RF_TED_mu_values_vec; RF_TED_mu_values];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conduct Comms Sampling
h = gen_firpm_first_order_h('low', 0.49*pi, 0.51*pi, 1, 30);
h = [h zeros(1, mod(-length(h), 2))];
ROM = genPolyphaseROM(h, 2, 'counter-clockwise');
[recieved_sampled_I, recieved_sampled_Q]= applyBarkerDecimator(ROM',Comms_I_symbols_recovered_no_farrow,Comms_Q_symbols_recovered_no_farrow,2);


recieved_sampled_symbol_axis_I = 0:1/Fs_in:length(recieved_sampled_I)/Fs_in - (1/Fs_in);
recieved_sampled_symbol_axis_Q = 0:1/Fs_in:length(recieved_sampled_Q)/Fs_in - (1/Fs_in);

% implement a rejection ratio to reject transients 
recieved_sampled_I = recieved_sampled_I(abs(recieved_sampled_I) >= rejection_ratio);
recieved_sampled_Q = recieved_sampled_Q(abs(recieved_sampled_Q) >= rejection_ratio);

% map the resulting symbols to data bits
try
    recovered_symbols = recieved_sampled_I + 1j*recieved_sampled_Q;
    recovered_bits = decMpsk(M, recovered_symbols, symbolMap);
catch
    recovered_bits = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consolidatae the output signals
Final_Recovered_Signals = {recovered_bits, RF_symbols_recovered_vec(1,:), RF_symbols_recovered_vec(2,:), RF_symbols_recovered_vec(3,:), RF_symbols_recovered_vec(4,:), RF_symbols_recovered_vec(5,:), RF_symbols_recovered_vec(6,:), RF_symbols_recovered_vec(7,:)};

for channel = 1:8
    fig_n = stemplot(fig_n, 1:length(Final_Recovered_Signals{channel}), Final_Recovered_Signals{channel}, ['Output Signal x' num2str(channel-1)], 'Index', 'Magnitude');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BER = sum(x0 ~= recovered_bits) / length(x0);
fprintf("Bit Error Rate: %f\n", BER);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%