% Sampling Practice
clear; clc; close all;


M = 16; % QAM order
numSymbols = 10000;
sps = 2; % Samples per symbol (oversampling)
rolloff = 0.35;
span = 8; % Filter span in symbols

% Generate random QAM symbols
data = randi([0 M-1], numSymbols, 1);
txSymbols = qammod(data, M, 'UnitAveragePower', true);

% Pulse shaping filter
rrcFilter = rcosdesign(rolloff, span, sps, 'sqrt');
txWaveform = upfirdn(txSymbols, rrcFilter, sps, 1);

% Add timing offset (fractional delay)
timingOffset = 0.3; % In symbol periods
rxWaveform = interp1(1:length(txWaveform), txWaveform, ...
    (1:length(txWaveform)) + timingOffset * sps, 'linear', 0);

% Add AWGN
snr = 20; % dB
rxWaveform = awgn(rxWaveform, snr, 'measured');

rxFiltered = filter(rrcFilter, 1, rxWaveform);












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = gardner_ted(r)
% r: complex baseband samples at 2 samples/symbol
% e: timing error signal (length: numSymbols-2)

    late = r(3:2:end);
    early = r(1:2:end-2);
    prompt = r(2:2:end-1);

    e = real((late - early) .* conj(prompt));
end


function [mu, state] = pi_loop_filter(e, state, K1, K2)
% e: timing error input
% state: struct with .integrator
% K1: proportional gain
% K2: integrator gain
% mu: updated timing phase

    state.integrator = state.integrator + K2 * e;
    mu = K1 * e + state.integrator;
end
