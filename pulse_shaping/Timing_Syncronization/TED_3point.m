% Root Cosine Pulse Shaping
close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare system parameters
L = 10;            % upsamples per symbol
M = 2;             % M-PSK
span = 10;          % Filter span (pulse shape filter order is span*L)
rolloff = 1;       % Alpha (0 goes to sinc, 1 foes to more square shaped in time) "excess bandwidth"
numSymbols = 100;  % Number of symbols in transmission 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate MPSK symbols, filter, upsample and filter
[bitGroups, symbolMap, MpskSymbols] = genMpsk(M, numSymbols);
h = genPulseFilter(rolloff, span, L, "sqrt");
upsampledSignal = upfirdn(MpskSymbols, 1, L);
upsampledFilterdSignal = filter(h, 1, upsampledSignal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot your filtered signals
figure;
subplot(2,1,1);
plot(upsampledSignal, 'LineWidth', 1.2);
title('Upsampled Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(upsampledFilterdSignal, 'LineWidth', 1.2);
title('Upsampled-Filtered Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

figure;
eyediagram(upsampledFilterdSignal, L);
title('Eye Diagram of Upsampled-Filtered Signal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% CHANNEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Gauss noise vector
gaussNoise = randn(length(upsampledFilterdSignal),1)'/5;
recievedSignal = upsampledFilterdSignal + gaussNoise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the Signal by 2 for Gardner TED
recievedInterpolatedSignal = upfirdn(recievedSignal, 1, 2);
% Apply the matched filter
matchFilteredSignal = filter(h, 1, recievedInterpolatedSignal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLot the recieved and post matched filter signals
figure;
subplot(2,1,1);
plot(recievedInterpolatedSignal, 'LineWidth', 1.2);
title('Received Interpolated Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(matchFilteredSignal, 'LineWidth', 1.2);
title('Matched-Filtered Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

figure;
eyediagram(matchFilteredSignal, L*2);
title('Eye Diagram of Matched Filtered Signal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Gardner Timing Error Detector
TED_Shift_reg = zeros(1,3);
lengthMatchFilteredSignal = length(matchFilteredSignal);
error = zeros(1, ceil(lengthMatchFilteredSignal/2));
phaseAdjustment = zeros(1, ceil(lengthMatchFilteredSignal/2));

for k = 1:lengthMatchFilteredSignal
    % Update the shift register with the received signal
    TED_Shift_reg = [matchFilteredSignal(k), TED_Shift_reg(1:2)];
    % Calculate the timing error
    timingError = (TED_Shift_reg(1) - TED_Shift_reg(3)) * TED_Shift_reg(2);
    phaseError = atan((TED_Shift_reg(1) - TED_Shift_reg(3))/(2*L*2));
    % Decimate the error calculation by 2 to get 1 error for each symbol
    if(mod(k,2) == 0)
        error(k/2) = timingError;
        phaseAdjustment(k/2) = phaseError;
    end
end

%Upsample Error signal to match the Interpolation clock
interpolatedError = upfirdn(error, 1, 2);
interpolatedPhaseAdjustment = upfirdn(phaseAdjustment, 1, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
stem(interpolatedError);
title('Gardner Timing Error Signal');
xlabel('Symbol Index');
ylabel('Error');
grid on;
figure;
stem(interpolatedPhaseAdjustment);
title('Gardner Timing Phase Signal');
xlabel('Symbol Index');
ylabel('Error rad/s');
grid on;













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Defs
function [bitGroups, symbolMap, MpskSymbols] = genMpsk(M, numSymbols)
    bits = randi([0 1], log2(M)*numSymbols, 1);
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
    h = h / max(h);
end