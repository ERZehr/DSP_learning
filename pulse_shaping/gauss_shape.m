% Gaussian Pulse Shaping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a random BPSK sequence
bpskSequence = randi([0 1], 1, 100); % Generate a random BPSK sequence
bpskSignal = 2*bpskSequence - 1; % Map 0 to -1 and 1 to 1
figure(5);
stem(bpskSignal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare and calculate Parameters
N = 2; % number of Int CIC stages
R = 20; % Interpolation Rate
M = 1; % Differential Delay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate impulse response coefficients
h = calc_impulse_response(N, R, M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3 : calculate the frequency response
[H, w] = freqz(h, 1, 4096);
H_dB = 20 * log10(abs(H));
H_dB = [flip(H_dB)' H_dB']';
w = [-flip(w)' w']';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the sequence by 10
interpolated_bpskSIgnal = upfirdn(bpskSignal, 1, 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply filter
PS_output = filter(h, 1, interpolated_bpskSIgnal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4 : Plot results
% Impulse Response M=1
figure(9);
stem(h, '-');
title('Gauss Impulse Response'); xlabel('Index'); ylabel('Magnitude'); grid on;
% Magnatude Response
figure(10);
plot(w, H_dB, '-');
title('Gauss Filter Frequency Response');
set(gca, 'XTick', -pi:pi/10:pi);
set(gca, 'XTickLabel', {'-pi', '-0.9pi', '-0.8pi', '-0.7pi', '-0.6pi', ...
    '-0.5pi', '-0.4pi', '-0.3pi', '-0.2pi', '-0.1pi', '0', '0.1pi', ...
    '0.2pi', '0.3pi', '0.4pi', '0.5pi', '0.6pi', '0.7pi', '0.8pi', ...
    '0.9pi', 'pi'});
xlabel('Normalized Frequency');
ylabel('Magnitude Response (dB)');
grid on;
% Filtered Data
figure(11);
plot(PS_output, '-');
hold on;
plot(PS_output_norm, '-');
title('Gauss Filter Pulse');
yline(1, '--'); yline(-1, '--');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Function
function h = calc_impulse_response(N, R, M)
    RM = R * M;
    base = ones(1, RM);
    h = base;
    for k = 2:N
        h = conv(h, base);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%