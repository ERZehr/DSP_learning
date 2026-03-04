% polyphase channelizer development
close all; clear; clc; fig_n = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channelizer Parameters
Fs_in = 100000; % 0-24kHz range
K = 8; % number of channels
M = 8; % number of polyphase components per channel
D = 8; % number of decimations
% Prototype Filter
wp = pi/K-pi/150;
ws = pi/K+pi/150;
rp = 0.1;
rs = 100;
filter_width = wp*Fs_in/pi;
center_frequency = filter_width/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Signal Declaration
t = 0:1/Fs_in:1-1/Fs_in;
L = length(t);

f_vec = [center_frequency center_frequency+filter_width center_frequency+(2*filter_width) center_frequency+(3*filter_width)]';
A_vec = [1 1 1 1]';
x1 = sum(A_vec.*cos(2*pi*f_vec.*t));
x2 = sum(A_vec.*sin(2*pi*f_vec.*t));

symbol_axis = 0:1/Fs_in:length(x2)/Fs_in - (1/Fs_in);
fig_n = contplot(fig_n, symbol_axis, x1, 'x1', 'Time (s)', 'Magnitude');
fig_n = contplot(fig_n, symbol_axis, x2, 'x2', 'Time (s)', 'Magnitude');

DFT1 = fft(x1);
DFT2 = fft(x2);
DFT1 = fftshift(DFT1);
DFT2 = fftshift(DFT2);

f_axis = (-L/2:L/2-1)*(Fs_in/L); % Hz

fig_n = fig_n + 1;
figure(fig_n);
plot(f_axis, abs(DFT1), 'LineWidth', 2)
title('Magnitude Spectrum of x1')
xlabel('Frequency (Hz)')
ylabel('|X1(f)|')
grid on

fig_n = fig_n + 1;
figure(fig_n);
plot(f_axis, abs(DFT2), 'LineWidth', 2)
title('Magnitude Spectrum of x2')
xlabel('Frequency (Hz)')
ylabel('|X2(f)|')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate prototype filter
% Plot prototype filter
h0 = gen_firpm_first_order_h('low', wp, ws, rp, rs);
h0 = [h0 zeros(1, mod(-length(h0), K))];
[H0,w] = freqz(h0,1,4096);         % Frequency response

% Generate channel filters
n = 0:length(h0)-1;
k = (0:K-1)';
modulated_vector = exp(1j*2*pi*(k+0.5)*n/K); % no negative here to match idft
shifted_filters = modulated_vector .* h0;

fig_n = fig_n + 1;
figure(fig_n);
hold on
for channel_index = 1:K
    H = freqz(shifted_filters(channel_index,:),1,4096,'whole');
    H = fftshift(H);
    w = linspace(-pi,pi,length(H));
    plot(w*Fs_in/2/pi,20*log10(abs(H)))
end
legend('1', '2', '3','4', '5', '6', '7', '8');
grid on
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
title('Channelized Filter Responses')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply shifted filters to inputs and plot
channel_outputs = [];
for k = 1:K
    channel_outputs(k,:) = filter(shifted_filters(k,:),1,x1);
end

for k = 1:K
    fig_n = fig_n + 1;
    figure(fig_n);
    plot(real(channel_outputs(k,:)))
    title('Channel %d Output', k)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice prototype filter into polyphase sections
% Plot polyphase sections
h0_poly = zeros(M, length(h0)/M);
for index = 1:M
    h0_poly(index,:) = h0(index:M:end);  % pick every K-th sample
end

fig_n = fig_n + 1;
figure(fig_n);
hold on;
for k = 1:size(h0_poly,1)
    taps = h0_poly(k,:);          % Extract row filter
    [H,f] = freqz(taps,1,1024);         % Frequency response
    plot(f/pi,20*log10(abs(H)));        % Plot magnitude in dB
end
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');
title('Frequency Responses of Polyphase Sections');
grid on;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate IDFT
lambda = 0:M-1;
k = (0:K-1)';
modulated_vector = exp(1j*2*pi*(k+0.5)*lambda/M); % no negative here to match idft

[M, L] = size(h0_poly);   % M polyphase branches, L taps per branch
channel_polyphase_coefficient_dot_products = zeros(K, M, L); % K channels × M branches × L taps
for index_k = 1:K
    for index_lambda = 1:M
        channel_polyphase_coefficient_dot_products(index_k, index_lambda,:) = modulated_vector(index_k, index_lambda) .* h0_poly(index_lambda,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply shifted filters to inputs and plot
channel_outputs = zeros(K,Fs_in+1);
for k = 1:K
    channel_outputs(k,:) = applyBarkerDecimator(channel_polyphase_coefficient_dot_products(k),x1,x2,1);
end

for k = 1:K
    fig_n = fig_n + 1;
    figure(fig_n);
    plot(real(channel_outputs(k,:)))
    title('Channel %d Output', k)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%