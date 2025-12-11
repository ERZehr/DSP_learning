% Sampling Practice
clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Discrete Signal Declaration
m = linspace(0, pi/20, 10000000);
x = cos(2*pi*400*m); % 400Hz Cosine

fig_n = 0;
fig_n = fig_n + 1;
figure(fig_n);
plot(m, x);
xlabel('Sample ');
ylabel('Magnitude');
title('x_m');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample to get x'_t
fs = 1000;      % sample rate 1 kHz
Ts = 1/fs;
N = 1000;       % number of samples to take
n = 0:N-1;      % sample indices

for k = 1:n
    sampled_signal(n) = x(k*Ts) * dirac(t-n*Ts);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%