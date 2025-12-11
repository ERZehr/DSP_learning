% Parameters
N = 4; R = 5; M = 1; input_bits = 7; acc_bits = 17;
num_samples = 100;

% Test inputs
n = 0:num_samples-1;
omega0 = 2*pi/500;
x_cos = cos(omega0 * n);
x_step = ones(1, num_samples);

% Quantize input for fixed-point
xq_cos = round(x_cos * (2^(input_bits-1)-1));
xq_cos = max(min(xq_cos, 2^(input_bits-1)-1), -2^(input_bits-1));
xq_step = round(x_step * (2^(input_bits-1)-1));

% Floating-point simulation
[stage_outs_fp_cos] = cic_filter_fp(x_cos, N, R, M);
[stage_outs_fp_step] = cic_filter_fp(x_step, N, R, M);

% Fixed-point simulation
[stage_outs_fx_cos] = cic_filter_fx(xq_cos, N, R, M, acc_bits);
[stage_outs_fx_step] = cic_filter_fx(xq_step, N, R, M, acc_bits);

% Plotting
plot_stage_outputs(stage_outs_fp_cos, 'Floating-Point Cosine');
plot_stage_outputs(stage_outs_fx_cos, 'Fixed-Point Cosine');
plot_stage_outputs(stage_outs_fp_step, 'Floating-Point Step');
plot_stage_outputs(stage_outs_fx_step, 'Fixed-Point Step');
