% DTFT, IDTFT, and FFT practice
clear; clc; close all;

Fs = 1000;
t = 0:1/Fs:1;
N = length(t);

f0 = 200;
x = cos(2*pi*f0*t);
%x = (0.9).^t;
[X_w_DTFT, w_DTFT] = DTFT(x, N);

X_w_FFT = fft(x, N);
f_FFT = [(-N/2:-1)*Fs/N (0:N/2)*Fs/N];

% Plot both signals
figure(1);
plot(w_DTFT, abs(X_w_DTFT), f_FFT, abs(X_w_FFT), 'LineWidth', 1.5);
xlabel('Frequency ');
ylabel('Magnitude');
title('DTFT of a signal');
grid on;



function [X, w] = DTFT(x, Nfreq)
    if nargin < 2
        Nfreq = 1024;
    end
    x = x(:)';
    n = 0:length(x)-1;
    w = linspace(-pi, pi, Nfreq);
    X = zeros(1, Nfreq);

    for k = 1:Nfreq
        X(k) = sum( x .* exp(-1j * w(k) * n) );
    end
end