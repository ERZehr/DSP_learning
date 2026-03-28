function [I_output, Q_output] = applyBarkerDecimator(ROM,I_input,Q_input,M)

    % Operate on input signal
    row_size = size(ROM,1);              % Number of taps in a ROM row
    shift_register = zeros(2, row_size); % Declare 2xrow size shift regester
    output_signal = [];                  % Declare 2xoutput vector
    symbols = [I_input; Q_input];        % Compress I and Q into a 2xlength matrix
    for k = 1:(size(symbols,2) + row_size)
        % Shift in the current input
        if k <= size(symbols,2)
            shift_register = [symbols(:,k) shift_register(:,1:end-1)]; % If we have actual values to shift in
        else
            shift_register = [zeros(2,1) shift_register(:,1:end-1)]; % If we are out and need to trail off
        end
        % Only produce output every M samples
        if mod(k-1, M) == 0
        
            % Sum over polyphase branches
            temp_sum = zeros(2,1);
            for p = 1:M
                % Dot product of each row with ROM(:,p)
                temp_sum(1) = temp_sum(1) + sum(shift_register(1,:) .* ROM(:,p).'); % multiply shift reg with rom and sum
                temp_sum(2) = temp_sum(2) + sum(shift_register(2,:) .* ROM(:,p).'); % multiply shift reg with rom and sum
            end
            % Correct column append
            output_signal = [output_signal temp_sum]; % append value to output register
        end
    end
    I_output = output_signal(1,:);
    Q_output = output_signal(2,:);
end