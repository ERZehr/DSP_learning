function [I_symbols_recovered,Q_symbols_recovered, I_symbols_recovered_no_farrow, Q_symbols_recovered_no_farrow, TED_error_vector, TED_mu_values] = apply_timing_data_recovery(I_input,Q_input,ROM_matched,tx_interp,Farrow_Length)

    row_size = size(ROM_matched,1);
    matched_filter_shift_register = zeros(2,row_size);
    %Decimation
    matched_filter_decimation_factor = tx_interp/2;
    % Gardner TED
    TED_mu_values = [];
    TED_error_vector = [];
    TED_error = [0;0];
    % Farrow Values
    mu = [1;1];
    farrow_sr = zeros(2,Farrow_Length);
    Farrow_A = vandermonde(Farrow_Length);
    polynomial_vals = zeros(1,Farrow_Length);

    symbols = [I_input; Q_input];

    % Output
    symbols_recovered = [];
    symbols_recovered_no_farrow = [];

    for k = 1:(size(symbols,2) + row_size) % ensure the whole signal runs through a zero initialized regester and runs all the way out
        % 1. Shift in single value
        if k <= size(symbols,2) % if we still have data to shift in 
            matched_filter_shift_register = [symbols(:,k) matched_filter_shift_register(:,1:end-1)]; % If we have actual values to shift in
        else                        % if we are out of data and trail off with zeros
            matched_filter_shift_register = [zeros(2,1) matched_filter_shift_register(:,1:end-1)]; % If we are out and need to trail off
        end

        % 2. Calculate output values
        % Only produce output every matched_filter_decimation_factor samples
        if mod(k-1, matched_filter_decimation_factor) == 0 % 2.1 Apply matched filter and decimate
            temp_sum = zeros(2,1);
            for p = 1:matched_filter_decimation_factor
                % Dot product of each row with ROM(:,p)
                temp_sum(1) = temp_sum(1) + sum(matched_filter_shift_register(1,:) .* ROM_matched(:,p).'); % multiply shift reg with rom and sum
                temp_sum(2) = temp_sum(2) + sum(matched_filter_shift_register(2,:) .* ROM_matched(:,p).'); % multiply shift reg with rom and sum
            end
            farrow_sr(1,:) = [temp_sum(1) farrow_sr(1,1:Farrow_Length-1)]; % Shift calculated I values into farrow shift register
            farrow_sr(2,:) = [temp_sum(2) farrow_sr(2,1:Farrow_Length-1)]; % Shift calculated Q values into farrow shift register
            % Polynomial branches
            for m = 1:Farrow_Length
                polynomial_vals(1,m) = Farrow_A(m,:) * farrow_sr(1,1:end)';
                polynomial_vals(2,m) = Farrow_A(m,:) * farrow_sr(2,1:end)';
            end
            % Apply mu powers
            Farrow_out = zeros(2,1);
            for m = 1:Farrow_Length % Calculate farrow output IQ values
                Farrow_out(1) = Farrow_out(1) + (mu(1)^(m-1))*polynomial_vals(1,m); % I factors
                Farrow_out(2) = Farrow_out(2) + (mu(2)^(m-1))*polynomial_vals(2,m); % Q factors
            end
            symbols_recovered(:,1:end+1) = [symbols_recovered(:,1:end) Farrow_out]; % Shift out altered value
            symbols_recovered_no_farrow(:,1:end+1) = [symbols_recovered_no_farrow(:,1:end) temp_sum]; % keep tabs on outputs without farrow filter
            % 2.2 Conduct Gardner TED calculation
            if(size(symbols_recovered,2)>2) % We need at least 3 values to run Gardner TED
                TED_error = symbols_recovered(:,end-1).*(symbols_recovered(:,end) - symbols_recovered(:,end-2)); % Calculate and normalize TED error for I
                mu = 1 + TED_error.*0.0001; % Apply error
            end

            TED_error_vector = [TED_error_vector TED_error]; % Log for plotting
            TED_mu_values = [TED_mu_values mu]; % Update mu log for plotting
        end
    end
    
    I_symbols_recovered = symbols_recovered(1,:);
    Q_symbols_recovered = symbols_recovered(2,:);
    I_symbols_recovered_no_farrow = symbols_recovered_no_farrow(1,:);
    Q_symbols_recovered_no_farrow = symbols_recovered_no_farrow(2,:);
end