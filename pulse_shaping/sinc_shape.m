% Sinc Pulse Shaping
close all; clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a random BPSK sequence
bpskSequence = randi([0 1], 1, 100); % Generate a random BPSK sequence
bpskSignal = 2*bpskSequence - 1; % Map 0 to -1 and 1 to 1
figure(1);
stem(bpskSignal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a windowed sinc function
fs = 10; % Sampling frequency
t = -4:1/fs:4; % Time vector (implicit rectangular window)
h = sinc(t); % Generate sinc pulse

[H, w] = freqz(h, 1, 4096);
H_dB = 20 * log10(abs(H));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the sequence by 10
interpolated_bpskSIgnal = upfirdn(bpskSignal, 1, 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply filter
PS_output = filter(h, 1, interpolated_bpskSIgnal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % plot results
 % Impulse Response
figure(2);
stem(h);
title('Sinc Impulse Response'); xlabel('Index'); ylabel('Magnitude'); grid on;
% Magnatude Response
figure(3);
plot(w, H_dB, '-');
title('Sinc Filter Frequency Response');
set(gca, 'XTick', -pi:pi/10:pi);
set(gca, 'XTickLabel', {'0', '0.1pi', '0.2pi', '0.3pi', '0.4pi', '0.5pi', ...
    '0.6pi', '0.7pi', '0.8pi', '0.9pi', 'pi'});
xlabel('Normalized Frequency');
ylabel('Magnitude Response (dB)');
grid on;
% Filtered Data
figure(4);
plot(PS_output); hold on;
title("Sinc Pulse Output Signals");
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%