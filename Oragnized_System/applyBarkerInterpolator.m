function [I_output,Q_output] = applyBarkerInterpolator(ROM,I_input,Q_input,L,C)

    phi_int = C/L;    % NCO Step
    accumulator = 0;

    [row_size, col_size] = size(ROM);
    shift_register = zeros(2,row_size);
    output_signal = [];
    symbols = [I_input; Q_input]; 

    for k = 1:(size(symbols,2) + row_size) % ensure the whole signal runs through a zero initialized regester and runs all the way out
        if k <= size(symbols,2)
            shift_register = [symbols(:,k) shift_register(:,1:end-1)];
        else
            shift_register = [zeros(2,1) shift_register(:,1:end-1)];
        end
        temp_sum = zeros(2,1);
        while accumulator < C
            ROM_addr_p = floor(col_size*accumulator/C);
            temp_sum(1) = sum(shift_register(1,:) .* ROM(:,ROM_addr_p+1).');
            temp_sum(2) = sum(shift_register(2,:) .* ROM(:,ROM_addr_p+1).');
            output_signal = [output_signal temp_sum];
            accumulator = accumulator + phi_int;
        end
        accumulator = accumulator - C;
    end
    I_output = output_signal(1,:);
    Q_output = output_signal(2,:);
end