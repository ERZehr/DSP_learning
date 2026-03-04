function [I_output, Q_output] = polyphaseIDFTChannelizer(ROM,I_input,Q_input,M)

    % Operate on input signal
    [row_size col_size] = size(ROM);
    shift_register = zeros(2, col_size); % Declare 2xrow size shift regester
    output_signals = [];                  % Declare 2xrow_sizexoutput_length vector
    symbols = [I_input; Q_input];              % Compress I and Q into a 2xlength matrix
    for k = 1:(size(symbols,2) + col_size)
        % Shift in the current input
        if k <= size(symbols,2)
            shift_register = [symbols(:,k) shift_register(:,1:end-1)]; % If we have actual values to shift in
        else
            shift_register = [zeros(2,1) shift_register(:,1:end-1)]; % If we are out and need to trail off
        end
        % Only produce output every M samples
        if mod(k-1, M) == 0
        
            % Sum over polyphase branches
            this_calc = zeros(row_size,2);
            for p = 1:row_size
                % Dot product of each row with ROM(:,p)
                this_calc(p,1) = ifft(sum(shift_register(1,:) * ROM.')); % multiply shift reg with rom and take ifft
                this_calc(p,2) = ifft(sum(shift_register(2,:) * ROM.')); % multiply shift reg with rom and take ifft
            end
            % Correct column append
            output_signals = [output_signals this_calc]; % append value to output register
        end
    end
    I_output = squeeze(output_signals(1,:,:));
    Q_output = squeeze(output_signals(2,:,:));
end