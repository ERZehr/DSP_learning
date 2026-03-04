function [I_mix,Q_mix] = IQ_mixer(I,Q,carrier_frequency,IQ_sample_frequency)
    n = 0:length(I)-1;

    % Define I and Q
    cos_lo = cos(2*pi*carrier_frequency/IQ_sample_frequency*n);
    sin_lo = sin(2*pi*carrier_frequency/IQ_sample_frequency*n);

    % geneteate I and Q mixed Signals, select I
    I_mix = I .* cos_lo - Q .* sin_lo; % The actual real signal
    Q_mix = I .* sin_lo + Q .* cos_lo;
end