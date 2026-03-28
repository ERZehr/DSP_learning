function [I_output, Q_output] = applyPolyphaseFFTChannelizer(ROM, I_input, Q_input, P)
%
%  ROM  : M x L matrix — row p holds the coefficients of polyphase filter p
%  P    : M, number of channels (commutator modulus)
    M   = P;
    row_size = size(ROM, 1);   % must equal M
    col_size = size(ROM, 2);   % taps per polyphase branch (L)
    shift_register_in_I = zeros(row_size, col_size);  % delay lines, one per branch
    shift_register_in_Q = zeros(row_size, col_size);  % delay lines, one per branch
    polyOut_I           = zeros(row_size, 1);          % one scalar per branch
    polyOut_Q           = zeros(row_size, 1);          % one scalar per branch
    shift_register_out_I = [];
    shift_register_out_Q = [];


    for k = 1:length(I_input) % Assuming I and Q have same length
        % --- Step 1: Commutator routes sample to branch p ---
        p = mod(k-1, M) + 1;   % 1-based branch index
        shift_register_in_I(p, :) = [I_input(k), shift_register_in_I(p, 1:end-1)];
        shift_register_in_Q(p, :) = [Q_input(k), shift_register_in_Q(p, 1:end-1)];
        % --- Step 2: Polyphase filter for branch p ---
        % convolve rom with input registers
        polyOut_I(p) = sum(ROM(p, :) .* shift_register_in_I(p, :));
        polyOut_Q(p) = sum(ROM(p, :) .* shift_register_in_Q(p, :));
        % --- Step 3: FFT across all M branches, once per cycle ---
        if p == M
            fftOut_I = fft(polyOut_I);               % M-point DFT → all channels
            fftOut_Q = fft(polyOut_Q);               % M-point DFT → all channels
            shift_register_out_I = [shift_register_out_I; fftOut_I.'];
            shift_register_out_Q = [shift_register_out_Q; fftOut_Q.'];
        end
    end
    % Each row of shift_register_out is one time instant; each column is a channel
    I_output = real(shift_register_out_I)';
    Q_output = real(shift_register_out_Q)';
end