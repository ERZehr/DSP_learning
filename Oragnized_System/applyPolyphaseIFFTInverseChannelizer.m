function [I_output, Q_output] = applyPolyphaseIFFTInverseChannelizer(ROM, I_input, Q_input, P)
%
%  ROM  : M x L matrix — row p holds the coefficients of polyphase filter p
%  P    : M, number of channels (commutator modulus)
    row_size = size(ROM, 1);   % must equal M
    col_size = size(ROM, 2);   % taps per polyphase branch (L)
    polyOut_I           = zeros(row_size, 1);          % one scalar per branch
    polyOut_Q           = zeros(row_size, 1);          % one scalar per branch
    shift_register_in_I = zeros(row_size, col_size);
    shift_register_in_Q = zeros(row_size, col_size);
    shift_register_out_I = [];
    shift_register_out_Q = [];


    for k = 1:size(I_input, 2) % Assuming I and Q have same length
        ifftOut_I = ifft(I_input(:,k));  % P-point DFT → all channels
        ifftOut_Q = ifft(Q_input(:,k));  % P-point DFT → all channels

        for s = 1:P
            shift_register_in_I(s,:) = [ifftOut_I(s) shift_register_in_I(s, 1:end-1)];
            shift_register_in_Q(s,:) = [ifftOut_Q(s) shift_register_in_Q(s, 1:end-1)];

            polyOut_I(s) = sum(ROM(s, :) .* shift_register_in_I(s));
            polyOut_Q(s) = sum(ROM(s, :) .* shift_register_in_Q(s));

            shift_register_out_I = [shift_register_out_I, polyOut_I(s)];
            shift_register_out_Q = [shift_register_out_Q, polyOut_Q(s)];
        end
    end
    % Each row of shift_register_out is one time instant; each column is a channel
    I_output = real(shift_register_out_I);
    Q_output = real(shift_register_out_Q);
end