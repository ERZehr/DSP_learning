% Root Cosine Pulse Shaping
close all; clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a random BPSK sequence
bpskSequence = randi([0 1], 1, 100); % Generate a random BPSK sequence
bpskSignal = 2*bpskSequence - 1; % Map 0 to -1 and 1 to 1
figure(5);
stem(bpskSignal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a windowed sinc function
sps = 10; % Sampling frequency
span = 8;
rolloff = 0.8; %(alpha)
h       = rcosdesign(rolloff, span, sps, "normal");
h_norm  = h / max(abs(h));
h_sqrt      = rcosdesign(rolloff, span, sps, "sqrt");
h_sqrt_norm = h_sqrt / max(abs(h_sqrt));

[H, w] = freqz(h, 1, 4096);
H_dB = 20 * log10(abs(H));

[H_norm, w_norm] = freqz(h_norm, 1, 4096);
H_dB_norm = 20 * log10(abs(H_norm));

[H_sqrt, w_sqrt] = freqz(h_sqrt, 1, 4096);
H_sqrt_dB = 20 * log10(abs(H_sqrt));

[H_sqrt_norm, w_sqrt_norm] = freqz(h_sqrt_norm, 1, 4096);
H_dB_sqrt_norm = 20 * log10(abs(H_sqrt_norm));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the sequence by 10
interpolated_bpskSIgnal = upfirdn(bpskSignal, 1, 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply filter
PS_output = filter(h, 1, interpolated_bpskSIgnal);
PS_output_norm = filter(h_norm, 1, interpolated_bpskSIgnal);
PS_sqrt_output = filter(h_sqrt, 1, interpolated_bpskSIgnal);
PS_sqrt_output_norm = filter(h_sqrt_norm, 1, interpolated_bpskSIgnal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % plot results
 % Impulse Response
figure(6);
stem(h);
hold on;
stem(h_norm);
hold on;
stem(h_sqrt);
hold on;
stem(h_sqrt_norm);
legend('raised cosine', 'norm raised cosine', 'root raised cosine', 'norm root raised cosine','Location', 'northeast');
title('Root Cosine Impulse Response'); xlabel('Index'); ylabel('Magnitude'); grid on;
hold off;
% Magnatude Response
figure(7);
plot(w, H_dB, '-', w_norm, H_dB_norm, '-', w_sqrt, H_dB_sqrt_norm, '-', w_sqrt_norm, H_dB_sqrt_norm, '-');
title('Root Cosine Filter Frequency Response');
set(gca, 'XTick', -pi:pi/10:pi);
set(gca, 'XTickLabel', {'0', '0.1pi','0.2pi', '0.3pi', '0.4pi', '0.5pi', ...
    '0.6pi', '0.7pi', '0.8pi', '0.9pi', 'pi'});
xlabel('Normalized Frequency');
ylabel('Magnitude Response (dB)');
legend('raised cosine', 'norm raised cosine', 'root raised cosine', 'root norm raised cosine','Location', 'northeast');
grid on;
% Filtered Data
figure;
plot(PS_output); hold on;
plot(PS_output_norm);
plot(PS_sqrt_output, '--');
plot(PS_sqrt_output_norm, '--');
yline(1,'--'); yline(-1,'--');
title("Pulse-Shaped Output Signals");
legend("RC","Norm RC","RRC","Norm RRC","Location","best");
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%