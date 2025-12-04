% Sampling Practice
clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Input Signal Declaration
Fs = 1000;
t = 0:1/Fs:1;

f0 = 20;
x = cos(2*pi*f0*t);
n = 1:length(x);

fig_n = 0;
fig_n = fig_n + 1;
figure(fig_n);
stem(n, x);
xlabel('Sample ');
ylabel('Magnitude');
title('Input Signal');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation Specs and LPF Generation and Plotting
M = 3;
wc = pi/M;
sig1 = 0.01; sig2 = 0.001;
wp = 0.85*wc; ws = 1.15*wc;
beta1 = sig1/2; beta2 = sig2;
f = [wp/pi ws/pi];
a = [1 0];
dev = [beta1 beta2];
[N,fo,ao,w] = firpmord(f,a,dev);
N = N + mod(M - mod(N,M), M) - 1; % this is so the number of taps is always divisable by L, so our ROM matrix is rectangular without zero padding
h = firpm(N,fo,ao,w);
[Hz, w] = freqz(h,1, 8192);
Hz_dB = 20 * log10(abs(Hz));

Rp_dB = 20*log10(1-sig1);
Rs_dB = 20*log10(sig2);
fig_n = fig_n + 1;
figure(fig_n);
plot(w, Hz_dB, '-');
title('Lowpass Decimation FIR Filter Frequency Response');
xline(wp, ':r', sprintf('wpH(z) = %.4f', wp)); xline(ws, ...
':r', sprintf('wsH(z) = %.4f', ws));
yline(Rp_dB/2, ':r', sprintf('Rp = %.4f dB', Rp_dB/2));
yline(-Rp_dB/2, ':r', sprintf('Rp = %.4f dB', -Rp_dB/2));
yline(Rs_dB, ':r', sprintf('Rs = %.4f dB', Rs_dB));
set(gca, 'XTick', 0:pi/10:pi);
set(gca, 'XTickLabel', {'0', '0.1pi', '0.2pi', '0.3pi', '0.4pi', ...
'0.5pi', '0.6pi', '0.7pi', '0.8pi', '0.9pi', 'pi'});
xlabel('Normalized Frequency');
ylabel('Magnitude Response (dB)');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter input signal and plot
x_filtered = conv(h, x);
y1 = downsample(x_filtered, M);
ny1 = 1:length(y1);
fig_n = fig_n + 1;
figure(fig_n);
stem(ny1, y1);
xlabel('Index');
ylabel('Magnitude');
title('Convolution of Filter and Input Signal');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split into Polyphase Components
polyphase_coeffs = cell(1, M);
for k = 1:M
   polyphase_coeffs{k} = h(k:M:end);
end

fig_n = fig_n + 1;
figure(fig_n);
for k = 1:M
   subplot(M, 1, k);
   stem(polyphase_coeffs{k}, 'filled');
   title(['Polyphase Component E_{', num2str(k-1), '}']);
   ylabel('h[n]');
   grid on;
end
xlabel('Tap Index');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populate ROM
row_size = length(h)/M;
column_size = M;

ROM = zeros(row_size, column_size);
for k = 1:column_size
    ROM(1:length(polyphase_coeffs{k}), k) = polyphase_coeffs{k}(:);
end
disp(ROM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Operate on input signal
shift_register = zeros(row_size,1);
output_signal = zeros(1, ceil(length(x)/M));
out_idx = 1;

for k = 1:(length(x)+row_size-1)
    % Shift in the current input
    if k <= length(x)
        shift_register = [x(k); shift_register(1:end-1)];
    else
        shift_register = [0; shift_register(1:end-1)];
    end
    % Only produce output every M samples
    if mod(k-1, M) == 0
        % sum over polyphase branches
        temp_sum = 0;
        for p = 1:M
            temp_sum = temp_sum + sum(shift_register .* ROM(:,p));
        end
        output_signal(out_idx) = temp_sum;
        out_idx = out_idx + 1;
    end
end
y2 = output_signal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the output signal
ny2 = 1:length(y2);
fig_n = fig_n + 1;
figure(fig_n);
stem(ny2, y2);
xlabel('Index');
ylabel('Magnitude');
title('ROM Based Polyphase Interpolation');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%