function data = concatenate_data(data)
    % fill the remaining empty timesteps with the last
    % voltage value we were able to simulate; this should have the
    % effect of allowing us to continue, just returning a high cost
    first_zero_idx = find(~data,1,'first');
    last_nonzero_idx = first_zero_idx - 1;
    data(first_zero_idx:end) = data(last_nonzero_idx);
end