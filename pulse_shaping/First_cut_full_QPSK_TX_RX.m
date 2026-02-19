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

% DAC Parameters
DAC_Amplification = 100;

% Channel Parameters
linear_channel_attenuation = 1;   % greater number = greater attenuation
gauss_noise_level = 1;            % greater number = greater noise
salt_and_pepper_noise_level = 0;  % greater number = greater noise

% ADC Parameters
ADC_interpolation_factor = 8;

% ADC Decimator Parameters
ADC_decimator_factor = ADC_interpolation_factor;

% Reciever LO/Filter Parameters
LO_recovery_cutoff = carrier_frequency * ADC_interpolation_factor / ADC_decimator_factor;

% Sampler Parameters
rejection_ratio = 0.3;


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
% Generate and apply interpolating pulse shaping filter
h = genPulseFilter(rolloff, span, L, "sqrt");
[H, w] = freqz(h, 1, 4096);
H = [flip(H)' H']; w = [-flip(w)' w'];
H_dB = 20*log10(abs(H));

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

% Apply pulse shaping filter to I
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

% Apply pulse shaping filter to Q
shift_register = zeros(row_size,1);
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

% geneteate I and Q mixed Signals, select I
I_mix = I_symbols_shpaed .* cos_lo - Q_symbols_shpaed .* sin_lo; % The actual real signal
Q_mix = I_symbols_shpaed .* sin_lo + Q_symbols_shpaed .* cos_lo;

s_prime = I_mix;
fig_n = stemplot(fig_n, symbol_axis, s_prime, 'Transmitted Signal', 'Time (s)', 'Magnitude');

% For DAC visualization purposes "convert the signal to analog"
s_prime = DAC_Amplification*s_prime;
fig_n = contplot(fig_n, symbol_axis, s_prime, 'DAC Output', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%% CHANNEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Gauss noise vector
gaussNoise = randn(length(s_prime),1)'*gauss_noise_level;

% Generate Salt and Pepper noise vector
prob = 0.05; % probability of impulse
saltAndPepperNoise = zeros(1,length(s_prime));
impulses = rand(1,length(s_prime)) < prob;
saltAndPepperNoise(impulses) = (2*randi([0,1],1,sum(impulses))-1)'*salt_and_pepper_noise_level;

total_noise = gaussNoise + saltAndPepperNoise;

s_prime_recieved = s_prime / linear_channel_attenuation + total_noise;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Potential ADC Interpolation
h_adc_model = fir1(128, 1/ADC_interpolation_factor);
h_adc_model = h_adc_model/max(h_adc_model);

s_prime_recieved = upfirdn(s_prime_recieved, h_adc_model, ADC_interpolation_factor);

symbol_axis_adc = 0:1/Fs_out/ADC_interpolation_factor:length(s_prime_recieved)/ADC_interpolation_factor/Fs_out - (1/Fs_out/ADC_interpolation_factor);
fig_n = stemplot(fig_n, symbol_axis_adc, s_prime_recieved, 'ADC Output', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compensate for ADC Interpolator TODO: make this a shift register
% polyphase

decimated_adc_output = upfirdn(s_prime_recieved, h_adc_model, 1, ADC_decimator_factor);
decimated_symbol_axis = 0:1/Fs_out:length(decimated_adc_output)/Fs_out - (1/Fs_out);
fig_n = stemplot(fig_n, decimated_symbol_axis, decimated_adc_output, 'Compensated ADC Output', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split I and Q  TODO: add simulated shift from Gardner Feedback
len_cos = 0:length(decimated_adc_output)-1;
cos_lo_rec = cos(2*pi*carrier_frequency/Fs_out*len_cos);
sin_lo_rec = sin(2*pi*carrier_frequency/Fs_out*len_cos);

I_recovered = 2*(decimated_adc_output .* cos_lo_rec);
Q_recovered = -2*(decimated_adc_output .* sin_lo_rec);

% Filter at carrier frequency to eliminate undesirable trig identity components
% we really just have to filter out 2Fs, so pass anything less then 5 and
% cut at higher than 1. Thus, this can be a really simple PM filter
wp = 0.5*pi; ws = 0.95*pi;
Rp = 1; Rs = 30;
dev_pass = (10^(Rp/20)-1)/(10^(Rp/20)+1);
dev_stop = 10^(-Rs/20);
f = [wp/pi ws/pi];
a = [1 0];
dev = [dev_pass dev_stop];
[N,fo,ao,w] = firpmord(f,a,dev);
h_IQ_recovery = firpm(N,fo,ao,w);

I_recovered = filter(h_IQ_recovery, 1, I_recovered);
Q_recovered = filter(h_IQ_recovery, 1, Q_recovered);

fig_n = stemplot(fig_n, decimated_symbol_axis, I_recovered, 'Recieved I Channel Signal', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, decimated_symbol_axis, Q_recovered, 'Recieved Q Channel Signal', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply matched filter to I, interpolate by 2 for Gardner, Gardner TED, and
% Farrow filter

% Matched filter
ROM_matched = ROM;
matched_filter_shift_register = zeros(row_size,1);
%Interpolation
matched_filter_interpolation_factor = 2;
accumulator = 0;
phi_int_matched_interopolator = C/matched_filter_interpolation_factor;
% Gardner TED
TED_error = 0;
TED_vector_I = [0];
% mu control loop
TED_proportional_gain = 0.0000002;
TED_integrator_gain = 0.0000005;
TED_loop_integration = 0;
mu = 0;
TED_garnder_values_I = [0];
% Farrow Values
Farrow_shift_register_I = [];
Farrow_new_value = 0;
Farrow_Length = 5;
Farrow_A = calcFarrowA(Farrow_Length);
Farrow_powers = (0:Farrow_Length-1).';


% Output
I_symbols_recovered = [];

for k = 1:(length(I_recovered)+row_size-1) % ensure the whole signal runs through a zero initialized regester and runs all the way out
    % 1. Shift in single value
    if k <= length(I_recovered) % if we still have data to shift in 
        matched_filter_shift_register = [I_recovered(k); matched_filter_shift_register(1:end-1)];
    else                        % if we are out of data and trail off with zeros
        matched_filter_shift_register = [0; matched_filter_shift_register(1:end-1)];
    end

    % 2. Calculate output values
    while accumulator < C
        % 2.1 Apply matched filter and interpolate
        ROM_addr_p = floor(matched_filter_interpolation_factor*accumulator/C) + 1;                  % Which matched filter ROM row values to use
        output_sum = sum(matched_filter_shift_register .* ROM_matched(:,ROM_addr_p)); % Calculate a single output value
        I_symbols_recovered = [I_symbols_recovered, output_sum];                      % Shift out a single value  
        accumulator = accumulator + phi_int_matched_interopolator;                    % update the NCO
        % 2.2 Conduct Gardner TED calculation
        if(length(I_symbols_recovered)>2)                                           % We need at least 3 values to run Gardner TED
            TED_error = I_symbols_recovered(end-1)*(I_symbols_recovered(end) - ...
                I_symbols_recovered(end-2)) / sqrt(abs(I_symbols_recovered(end))^2 + ...
                abs(I_symbols_recovered(end-1))^2 + abs(I_symbols_recovered(end-1))^2);  % Calculate and normalize TED error
            TED_vector_I = [TED_vector_I TED_error];                                     % Log for plotting
            TED_loop_integration = TED_loop_integration + TED_integrator_gain*TED_error; % Apply integration gain
            mu = 1 - mu + TED_proportional_gain*TED_error + TED_loop_integration;            % Apply proportional gain
            TED_garnder_values_I = [TED_garnder_values_I mu];                            % Update mu log for plotting
        end

        % 2.3 Update Farrow Values
        if(length(I_symbols_recovered)>=Farrow_Length)       % We need at least N values to run a N point Farrow filter
            Farrow_new_value = sum((mu.^Farrow_powers) .* (Farrow_A * I_symbols_recovered((end-Farrow_Length+1):end)'));
            Farrow_shift_register_I = [Farrow_shift_register_I Farrow_new_value];
        else
            Farrow_shift_register_I = [Farrow_shift_register_I 0];
        end


    end
    accumulator = accumulator - C;
end
TED_vector_I = [TED_vector_I 0];
TED_garnder_values_I = [TED_garnder_values_I TED_garnder_values_I(end)];
gardner_rate_axis = 0:1/Fs_out/matched_filter_interpolation_factor:length(I_symbols_recovered)/matched_filter_interpolation_factor/Fs_out - (1/Fs_out/matched_filter_interpolation_factor);
fig_n = stemplot(fig_n, gardner_rate_axis, TED_vector_I, 'Gardner Values I', 'Time (s)', 'Magnitude');
fig_n = contplot(fig_n, gardner_rate_axis, TED_garnder_values_I, 'Gardner mu I', 'Time (s)', 'Magnitude', '--');
fig_n = stemplot(fig_n, gardner_rate_axis, Farrow_shift_register_I, 'Farrow Filter Output', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply matched filter to Q, interpolate by 2 for Gardner, Gardner TED, and
% Farrow filter

% Matched filter
matched_filter_shift_register = zeros(row_size,1);
%Interpolation
accumulator = 0;
% Gardner TED
TED_error = 0;
TED_vector_Q = [0];
% mu control loop
TED_loop_integration = 0;
mu = 0;
TED_garnder_values_Q = [0];
% Farrow Values
Farrow_shift_register_Q = [];
Farrow_new_value = 0;


% Output
Q_symbols_recovered = [];

for k = 1:(length(Q_recovered)+row_size-1) % ensure the whole signal runs through a zero initialized regester and runs all the way out
    % 1. Shift in single value
    if k <= length(Q_recovered) % if we still have data to shift in 
        matched_filter_shift_register = [Q_recovered(k); matched_filter_shift_register(1:end-1)];
    else                        % if we are out of data and trail off with zeros
        matched_filter_shift_register = [0; matched_filter_shift_register(1:end-1)];
    end

    % 2. Calculate output values
    while accumulator < C
        % 2.1 Apply matched filter and interpolate
        ROM_addr_p = floor(matched_filter_interpolation_factor*accumulator/C) + 1;                  % Which matched filter ROM row values to use
        output_sum = sum(matched_filter_shift_register .* ROM_matched(:,ROM_addr_p)); % Calculate a single output value
        Q_symbols_recovered = [Q_symbols_recovered, output_sum];                      % Shift out a single value  
        accumulator = accumulator + phi_int_matched_interopolator;                    % update the NCO
        % 2.2 Conduct Gardner TED calculation
        if(length(Q_symbols_recovered)>2)                                           % We need at least 3 values to run Gardner TED
            TED_error = Q_symbols_recovered(end-1)*(Q_symbols_recovered(end) - ...
                Q_symbols_recovered(end-2)) / sqrt(abs(Q_symbols_recovered(end))^2 + ...
                abs(Q_symbols_recovered(end-1))^2 + abs(Q_symbols_recovered(end-1))^2);  % Calculate and normalize TED error
            TED_vector_Q = [TED_vector_Q TED_error];                                     % Log for plotting
            TED_loop_integration = TED_loop_integration + TED_integrator_gain*TED_error; % Apply integration gain
            mu = 1 - mu + TED_proportional_gain*TED_error + TED_loop_integration;            % Apply proportional gain
            TED_garnder_values_Q = [TED_garnder_values_Q mu];                            % Update mu log for plotting
        end

        % 2.3 Update Farrow Values
        if(length(Q_symbols_recovered)>=Farrow_Length)       % We need at least N values to run a N point Farrow filter
            Farrow_new_value = sum((mu.^Farrow_powers) .* (Farrow_A * Q_symbols_recovered(end-Farrow_Length+1:end)'));
            Farrow_shift_register_Q = [Farrow_shift_register_Q Farrow_new_value];
        else
            Farrow_shift_register_Q = [Farrow_shift_register_Q 0];
        end

    end
    accumulator = accumulator - C;
end
TED_vector_Q = [TED_vector_Q 0];
TED_garnder_values_Q = [TED_garnder_values_Q TED_garnder_values_Q(end)];
fig_n = stemplot(fig_n, gardner_rate_axis, TED_vector_Q, 'Gardner Values Q', 'Time (s)', 'Magnitude');
fig_n = contplot(fig_n, gardner_rate_axis, TED_garnder_values_Q, 'Gardner mu Q', 'Time (s)', 'Magnitude', '--');
fig_n = stemplot(fig_n, gardner_rate_axis, Farrow_shift_register_Q, 'Farrow Filter Output', 'Time (s)', 'Magnitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample the recovered signal
recieved_sampled_I = Farrow_shift_register_I(1:matched_filter_interpolation_factor*F_int:end);
recieved_sampled_Q = Farrow_shift_register_Q(1:matched_filter_interpolation_factor*F_int:end);


recieved_sampled_symbol_axis_I = 0:1/Fs_in:length(recieved_sampled_I)/Fs_in - (1/Fs_in);
recieved_sampled_symbol_axis_Q = 0:1/Fs_in:length(recieved_sampled_Q)/Fs_in - (1/Fs_in);
fig_n = stemplot(fig_n, recieved_sampled_symbol_axis_I, recieved_sampled_I, 'Sampled I', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, recieved_sampled_symbol_axis_Q, recieved_sampled_Q, 'Sampled Q', 'Time (s)', 'Magnitude');

% implement a rejection ratio to reject transients 
threshold = max([recieved_sampled_I recieved_sampled_Q])*rejection_ratio;
recieved_sampled_I = recieved_sampled_I(abs(recieved_sampled_I) >= threshold);
recieved_sampled_Q = recieved_sampled_Q(abs(recieved_sampled_Q) >= threshold);

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

function firpm_h = firpm_first_order_filter_gen(pass_type, wp, ws, Rp, Rs)
    dev_pass = (10^(Rp/20)-1)/(10^(Rp/20)+1);
    dev_stop = 10^(-Rs/20);
    switch pass_type
        case 'low'
            f = [wp/pi ws/pi];
            a = [1 0];
            dev = [dev_pass dev_stop];
        case 'high'
            f = [ws/pi wp/pi];
            a = [0 1];
            dev = [dev_stop dev_pass];
        otherwise
            f = [wp/pi ws/pi];
            a = [1 0];
            dev = [dev_pass dev_stop];
    end
    [N,fo,ao,w] = firpmord(f,a,dev);
    firpm_h = firpm(N,fo,ao,w);
end

function bitsRecovered = decMpsk(M, MpskSymbols, symbolMap)
    % Pre-allocate recovered symbol indices
    numSymbols = length(MpskSymbols);
    recoveredIndices = zeros(numSymbols, 1);

    % Map each received symbol to the nearest constellation point
    for k = 1:numSymbols
        [~, idx] = min(abs(MpskSymbols(k) - symbolMap)); % nearest constellation point
        recoveredIndices(k) = idx; 
    end

    % Convert Gray code indices back to binary indices
    if M == 2
        % For BPSK, simply map +1 -> 1, -1 -> 0
        bitsRecovered = (real(MpskSymbols) > 0);
    else
        % Reverse Gray coding
        grayIndices = recoveredIndices - 1;          % MATLAB 1-index adjustment
        binIndices = zeros(size(grayIndices));
        for n = 1:length(grayIndices)
            g = grayIndices(n);
            b = 0;
            mask = g;
            while mask > 0
                b = bitxor(b, mask);
                mask = bitshift(mask, -1);
            end
            binIndices(n) = b;
        end

        % Convert integer symbol indices to bits
        bitsPerSymbol = log2(M);
        bitsRecovered = de2bi(binIndices, bitsPerSymbol, 'left-msb'); 
        bitsRecovered = reshape(bitsRecovered.', [], 1);  % column vector
    end
end


function Farrow_A = calcFarrowA(order)
    t = (-floor(order/2):floor(order/2))';
    V = zeros(order,order);
    for i=1:order
        for j=0:order-1
            V(i,j+1) = t(i)^j;
        end
    end
    Farrow_A = inv(V)';
end
