% Root Cosine Pulse Shaping
close all; clear; clc;

% Declare system parameters
L = 10;            % upsamples per symbol
M = 2;             % M-PSK
span = 10;          % Filter span (pulse shape filter order is span*L)
rolloff = 1;       % Alpha (0 goes to sinc, 1 foes to more square shaped in time)
numSymbols = 100;  % Number of symbols in transmission  


% Generate MPSK symbols
[bitGroups, symbolMap, MpskSymbols] = genMpsk(M, numSymbols);

%%
% Plot the MPSK constellation 
figure(1); % Mpsk constellation
scatter(real(symbolMap), imag(symbolMap), 'filled', 'MarkerFaceColor','b');
axis equal;
title('Random MPSK Constellation');
grid on;
hold off;


% Generate pulse shaping filter
h = genPulseFilter(rolloff, span, L, "normal");
[H, w] = freqz(h, 1, 4096);
H = [flip(H)' H']; w = [-flip(w)' w'];
H_dB = 20*log10(abs(H));


% Plot the impulse and frequency responses of the pulse shaping filter
figure(2); % Full Pulse Shaping filter impulse response
stem(h);
title('RRC Impulse Response');
xlabel('Index'); ylabel('Normalized Magnitude');
grid on;

figure(3); % Full Pulse Shaping filter frequency response
plot(w, H_dB);
title('RRC Frequency Response');
set(gca, 'XTick', -pi:pi/10:pi);
set(gca, 'XTickLabel', {'-pi', '-0.9pi', '-0.8pi', '-0.7pi', '-0.6pi', ...
    '-0.5pi', '-0.4pi', '-0.3pi', '-0.2pi', '-0.1pi', '0', '0.1pi','0.2pi', ...
    '0.3pi', '0.4pi', '0.5pi', '0.6pi', '0.7pi', '0.8pi', '0.9pi', 'pi'});
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
grid on;
 

% Break Pulse Shaping Filter into L polyphase components
polyphaseFilters = genPolyphase(L, h);


% Apply the polypohase filters to each unique symbol
uniqueShapedPulses = cell(1, length(symbolMap));
for k = 1:length(symbolMap)
    paddedUNiqueSymbolMap = [zeros(1,3*L) symbolMap(k) zeros(1,3*L)];
    uniqueShapedPulses{k} = applyPolyphaseFilters(paddedUNiqueSymbolMap, polyphaseFilters, L);
end


% Apply the polypohase filters to the whole signal
    dataVector = zeros(1, L)';
    filterSums = cell(1, L);
    filtershiftRegisters = cell(L, L);
    polyphaseFIltersRam = polyphaseFilters;
    outputRegister = zeros(1,length(MpskSymbols));
    indexReg = 0;

    for k = 1:length(MpskSymbols)
        dataVector(indexReg+1) = MpskSymbols(k);
        for p = 1:L
            for q = 1:L
                filterSums()
            end
        end
    end


[H_signal, w] = freqz(outputRegister, 1, 4096);
H_signal = [flip(H_signal)' H_signal']; w = [-flip(w)' w'];
H_signal_dB = 20*log10(abs(H_signal));

% Plot the summed output signal and its real and complex components
figure(6); % Summed
plot(outputRegister);
title("RRC Pulse Summed Output Signals");
yline(1,'--'); yline(-1,'--');
grid on;

figure(7); % Summed filter output frequency domain
plot(w, H_signal_dB);
title('Filter Output in Frequency Domain');
set(gca, 'XTick', -pi:pi/10:pi);
set(gca, 'XTickLabel', {'-pi', '-0.9pi', '-0.8pi', '-0.7pi', '-0.6pi', ...
    '-0.5pi', '-0.4pi', '-0.3pi', '-0.2pi', '-0.1pi', '0', '0.1pi','0.2pi', ...
    '0.3pi', '0.4pi', '0.5pi', '0.6pi', '0.7pi', '0.8pi', '0.9pi', 'pi'});
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
grid on;

figure(8); % Real
hold on;
plot(real(outputRegister));
title("RRC Summed Real Output Signals");
yline(1,'--'); yline(-1,'--');
grid on;
hold off;

figure(9);% Complex
plot(imag(outputRegister));
title("RRC Summed Complex Output Signals");
yline(1,'--'); yline(-1,'--');
grid on;
hold off;



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

function [outputSignal] = applyPolyphaseFilters(signal, polyphaseFilters, L)
    shiftReg = zeros(1,L);
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