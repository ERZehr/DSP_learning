% polyphase channelizer development with Andrew's book
close all; clear; clc; fig_n = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channelizer Parameters
Fs_in = 100000; % 0-24kHz range
K = 8; % number of channels
D = 8; % number of decimations
% Prototype Filter
wp = pi/K-pi/100;
ws = pi/K+pi/100;
rp = 0.1;
rs = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Signal Declaration
t = 0:1/Fs_in:1-1/Fs_in;
L = length(t);

f_vec = [100 1000 2500 10000]';
A_vec = [1 7 3 1]';
x1 = sum(A_vec.*cos(2*pi*f_vec.*t));
%x1 = cos(2*pi*100.*t);
x2 = sum(A_vec.*sin(2*pi*f_vec.*t));
%x2 = sin(2*pi*100.*t);

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
modulated_vector = exp(1j*2*pi*(k)*n/K); % no negative here to match idft
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
legend(string(1:K));
grid on
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
title('Channelized Filter Responses')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Polyphase FFT Channelizer
ROM = genPolyphaseROM(h0, K, 'clockwise');
[I_output, Q_output] = applyPolyphaseFFTChannelizer(ROM, x1, x2, K);

Ts = 1/(Fs_in * K);
symbol_axis = (0:size(I_output,2)-1) * Ts;

% Loop over channels and plot each one
for ch = 1:K/2+1
    fig_n = stemplot(fig_n, symbol_axis, I_output(ch,:), ['I\_output Channel ' num2str(ch)], 'Time (s)', 'Magnitude');
end
for ch = 1:K/2+1
    fig_n = stemplot(fig_n, symbol_axis, Q_output(ch,:), ['Q\_output Channel ' num2str(ch)], 'Time (s)', 'Magnitude');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Inverse Polyphase FFT Channelizer
ROM = genPolyphaseROM(h0, K, 'counter-clockwise');
[I_output, Q_output] = applyPolyphaseIFFTInverseChannelizer(ROM, I_output, Q_output, K);

I_output = K^2 * I_output;
Q_output = K^2 * Q_output;

Ts = 1/(Fs_in);
symbol_axis = (0:length(I_output)-1) * Ts;

% Plot output channels
fig_n = stemplot(fig_n, symbol_axis, x1, 'x1', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, I_output, 'I_output Channel', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, x2, 'x2', 'Time (s)', 'Magnitude');
fig_n = stemplot(fig_n, symbol_axis, Q_output, 'Q_output Channel', 'Time (s)', 'Magnitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%