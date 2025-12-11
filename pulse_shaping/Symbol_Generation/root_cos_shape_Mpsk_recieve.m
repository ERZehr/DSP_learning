% Root Cosine Pulse Shaping
close all; clear; clc;

% Declare system parameters
L = 10;            % upsamples per symbol
M = 4;             % M-PSK
span = 8;          % Filter span (pulse shape filter order is span*L)
rolloff = 0.2;       % Alpha (0 goes to sinc, 1 foes to more square shaped in time)
numSymbols = 1000;  % Number of symbols in transmission  

% Generate MPSK symbols
[bitGroups, symbolMap, MpskSymbols] = genMpsk(M, numSymbols);
% Generate pulse shaping filter
h = genPulseFilter(rolloff, span, L, "sqrt");
% Filter the full Mpsk sequence with the Pulse Shaping filter
upsampledSignal = upfirdn(MpskSymbols, 1, L);
upsampledFilterdSignal = filter(h, 1, upsampledSignal);


% Generate Gauss noise vector
noiseLength = 50000;
gaussNoise = randn(noiseLength,1)';
figure(1)
plot(gaussNoise);

%generate a random index in the gauss noise for my signal to "come back"
startDelayLength = 30000;
startIndex = randi(startDelayLength);


% generate the signal that comes back
% recievedSignal = gauss noise + signal at the random time that it comes back
delayedSignal = [zeros(1,startIndex) upsampledFilterdSignal zeros(1, noiseLength - (startIndex + length(upsampledFilterdSignal)))];
idx = find(delayedSignal ~= 0, 1, 'first'); % check where the signal starts
fprintf('Actual Delay Index = %d\n', idx);
% Plot the delayed signal
figure(2);
plot(delayedSignal);
title('Delayed Signal');
xlabel('Index'); ylabel('Magnitude');
grid on;


% Calculate the recieved data
% recievedSignal = gauss noise + signal at the random time that it comes back
recievedSignal = gaussNoise + delayedSignal;
% Plot the delayed function
figure(3);
plot(recievedSignal);
title('Recieved Signal');
xlabel('Index'); ylabel('Magnitude');
grid on;


% Calculate the autocorrelation of the sequence
autocorreleation = zeros(1,length(recievedSignal)-length(upsampledFilterdSignal) + 1);
for k = 1:length(autocorreleation)
    recievedChunk = recievedSignal(k:k+length(upsampledFilterdSignal)-1);
    autocorreleation(k) = sum(upsampledFilterdSignal .* recievedChunk);
end

% Find the max of autocorrelation
[autocorreleationMaxValue, autocorreleationMaxIndex]  = max(autocorreleation);
fprintf("Calculated Delay Index = %d", autocorreleationMaxIndex);

% Plot the autocorrelation function
figure(4);
plot(autocorreleation);
title('Autocorrelation of sent and recieved data');
xlabel('Index'); ylabel('Normalized Magnitude');
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