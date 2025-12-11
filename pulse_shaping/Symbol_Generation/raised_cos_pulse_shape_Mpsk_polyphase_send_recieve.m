% Root Cosine Pulse Shaping
close all; clear; clc;

% Declare system parameters
L = 10;            % upsamples per symbol
D = 10;            % decimations per sample
M = 2;             % M-PSK
span = 10;          % Filter span (pulse shape filter order is span*L)
rolloff = 1;       % Alpha (0 goes to sinc, 1 foes to more square shaped in time)
numSymbols = 100;  % Number of symbols in transmission  

% Generate MPSK symbols
[bitGroups, symbolMap, MpskSymbols] = genMpsk(M, numSymbols);
figure(1); % Full Pulse Shaping filter impulse response
stem(MpskSymbols);
title('Original Signal');
xlabel('Index'); ylabel('Normalized Magnitude');
grid on;

% Generate pulse shaping filter
h = genPulseFilter(rolloff, span, L, "normal");
figure(2); % Full Pulse Shaping filter impulse response
stem(h);
title('Impulse Response');
xlabel('Index'); ylabel('Normalized Magnitude');
grid on;

% Break Pulse Shaping Filter into L polyphase components
polyphaseFilters = genPolyphase(L, h);

% Apply the polypohase filters to the whole signal
filteredSignal = polyphaseInterpolator(MpskSymbols, polyphaseFilters, L);
figure(3); % Full Pulse Shaping filter impulse response
plot(filteredSignal);
title('Filtered Signal');
xlabel('Index'); ylabel('Normalized Magnitude');
grid on;

% Generate Channel Noise
channelNoise = randn(numSymbols*L,1)'/3;
channelSignal = filteredSignal + channelNoise;
figure(4); % Full Pulse Shaping filter impulse response
plot(channelSignal);
title('Channel Signal');
xlabel('Index'); ylabel('Normalized Magnitude');
grid on;

unfilteredSignal = polyphaseDecimator(channelSignal, polyphaseFilters, D);
figure(5); % Full Pulse Shaping filter impulse response
plot(unfilteredSignal);
title('Unfiltered Signal');
xlabel('Index'); ylabel('Normalized Magnitude');
grid on;

recoveredMpskBitGroups = decodeSymbols(unfilteredSignal, M);
figure(6); % Full Pulse Shaping filter impulse response
stem(recoveredMpskBitGroups);
title('Recovered Signal');
xlabel('Index'); ylabel('Normalized Magnitude');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function [outputSignal] = decodeSymbols(signal, M)
    theta = angle(signal);
    q = M/(2*pi) * theta;
    outputSignal = mod(round(q), M);
    if M==2
        outputSignal = -(2*outputSignal - 1)';
    else
        fprintf("a");
    end
end

function [h] = genPulseFilter(rolloff, span, L, type)
    h = rcosdesign(rolloff, span, L, type);
    h = h / max(h);
end

function [polyphaseFilters] = genPolyphase(L, h)
    polyphaseFilters = cell(1, L);
    for k = 1:L
        polyphaseFilters{k} = h(k:L:end);
    end
end

function [outputSignal] = polyphaseInterpolator(signal, polyphaseFilters, L)
    branchOutputs  = cell(1, length(polyphaseFilters));
    lengths = zeros(1, length(polyphaseFilters));
    % Apply Filters and upsample
    for k = 1:length(polyphaseFilters)
        branchOutputs{k} = filter(polyphaseFilters{k}, 1, signal);
        branchOutputs{k} = upfirdn(branchOutputs{k}, 1, L);
    end
    % Apply delay values for polyphase reconstruction
    for k = 1:length(polyphaseFilters)
        branchOutputs{k} = [zeros(1,k-1) branchOutputs{k}];
        lengths(k) = length(branchOutputs{k});
    end
    % Find the max length of branch outputs to pad if (taps%L != 0)
    maxLen = max(lengths);
    for k = 1:length(polyphaseFilters)
        branchOutputs{k}(end+1:maxLen) = 0;
    end
    % Combine the delayed polyphase signals back together
    outputSignal = sum(cat(1, branchOutputs{:}), 1);
end

% There is a bit of a bug in this one, I think it has to do with 
% transient behavior, I get D samples of garbage prior to the signal
% then the last D samples are a bit shoddy. 
function [outputSignal] = polyphaseDecimator(signal, polyphaseFilters, D)
    branchOutputs  = cell(1, length(polyphaseFilters));
    % Delay the input signals
    for k = 1:length(polyphaseFilters)
        branchOutputs{k} = [zeros(1, k) signal];
    end
    % Decimate the delayed signals
    for k = 1:length(polyphaseFilters)
        branchOutputs{k} = branchOutputs{k}(1:D:end);
    end
    %Filter the signals
    for k = 1:length(polyphaseFilters)
        branchOutputs{k} = filter(polyphaseFilters{k}, 1, branchOutputs{k});
    end
    % Combine the polyphase signals back together
    outputSignal = sum(cat(1, branchOutputs{:}), 1);
    outputSignal = outputSignal(D+1:end);
end