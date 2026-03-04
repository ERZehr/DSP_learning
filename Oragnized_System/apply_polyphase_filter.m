function y = apply_polyphase_filter(x, h_poly)

    % Preallocate cell array for branch outputs
    M = size(h_poly, 1);
    branch_outputs = cell(1, M);

    % Apply each polyphase branch
    for branch_idx = 1:M
        branch_outputs{branch_idx} = upfirdn(x, h_poly(branch_idx,:), 1, M);
    end

    % Sum all branch outputs
    y_length = max(cellfun(@length, branch_outputs));
    y = zeros(1, y_length);

    for branch_idx = 1:M
        % zero-pad branch output if necessary
        branch = branch_outputs{branch_idx};
        if length(branch) < y_length
            branch = [branch zeros(1, y_length - length(branch))];
        end
        y = y + branch;
    end
end