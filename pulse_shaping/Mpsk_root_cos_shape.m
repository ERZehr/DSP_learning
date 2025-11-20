% Root Cosine Pulse Shaping
close all; clear; clc;

L = 10;          % Samples per symbol
M = 8;           % M-PSK
span = 8;        % Filter span
rolloff = 0.5;     % Alpha
numSymbols = 1000;

% Generate M-PSK symbols
bits = randi([0 1], log2(M)*numSymbols, 1);
bitGroups = reshape(bits, log2(M), []).';
symbolMap = exp(1j * 2 * pi * (0:M-1) / M);
indices = bi2de(bitGroups, 'left-msb') + 1;
MpskSymbols = symbolMap(indices);

figure(1);
scatterplot(MpskSymbols);
title('Random MPSK Constellation');
grid on;
hold off;

% RC pulse shaping filter
h = rcosdesign(rolloff, span, L, "sqrt");

[H, w] = freqz(h, 1, 4096);
H_dB = 20*log10(abs(H));

figure(2);
stem(h);
title('RRC Impulse Response');
grid on;

figure(3);
plot(w, H_dB);
title('RRC Frequency Response');
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
grid on;

% Interpolation
upsampled = upfirdn(MpskSymbols, 1, L);

% Pulse shaping
PS_output = filter(h, 1, upsampled);

% Plot shaped waveform
figure(4);
plot(PS_output);
title("Pulse-Shaped Output Waveform");
grid on;

% Eye diagram
figure(5);
eyediagram(PS_output, L);
title('Eye Diagram (after Pulse Shaping)');