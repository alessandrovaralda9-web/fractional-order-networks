%% Simple Example: Fractional-Order Network Parameter Estimation
% This script demonstrates the basic usage of the fractional-order network
% parameter estimation code with a small example network.
%
% Author: [Your Name]
% Date: January 2025

clearvars;
close all;
clc;

fprintf('========================================\n');
fprintf('  Fractional-Order Network Example\n');
fprintf('========================================\n\n');

%% Step 1: Generate a Small Test Network
fprintf('Step 1: Generating a small test network...\n');

% Network dimensions (small for quick demonstration)
n = 6;   % Number of states
m = 2;   % Number of controllable inputs
v = 1;   % Number of uncontrollable inputs  
o = 2;   % Number of outputs

% Generate a single sparse network
[A, B, C, D, G, alpha, x0] = generateSparseNetwork(n, m, v, o, 'high', 1);

fprintf('  - States: %d\n', n);
fprintf('  - Controllable inputs: %d\n', m);
fprintf('  - Uncontrollable inputs: %d\n', v);
fprintf('  - Outputs: %d\n', o);
fprintf('  - Fractional orders (alpha): [%.2f, %.2f, %.2f, ...]\n', alpha(1), alpha(2), alpha(3));

%% Step 2: Create FON Object and Perform Decomposition
fprintf('\nStep 2: Creating FON object and decomposing...\n');

% Create state-space system and FON object
sys = ss(A, B, C, D);
fon = FON_Graph_Class(sys, alpha, G);

% Perform hierarchical decomposition
sub_systems = fon.hierarchicalDecomposition();
fprintf('  - Network decomposed into %d subsystems\n', length(sub_systems));

%% Step 3: Simulate the True System
fprintf('\nStep 3: Simulating the true system...\n');

% Simulation parameters
N = 50;   % Number of time steps (reduced for quick demo)
J = 50;   % Memory length

% Generate inputs
rng(42);  % For reproducibility
u = ones(m, N);  % Step inputs
w = 0.3 * randn(v, N);  % Random disturbances
x0_scaled = x0 * 100;  % Scale initial condition
e = 0.1 * randn(o, N);  % Small measurement noise

% Simulate
fon_true = FON_Graph_Class(sys, alpha, G);
fon_true.fsim(u, w, x0_scaled, e, J);

fprintf('  - Simulation completed with %d time steps\n', N);

%% Step 4: Parameter Estimation - Stage 1 (A and alpha)
fprintf('\nStep 4: Estimating A matrix and fractional orders...\n');

% Zero-input simulation for A and alpha estimation
u_zero = zeros(m, N);
w_zero = 0.3 * randn(v, N);
e_zero = 0.1 * randn(o, N);

fon_no_input = FON_Graph_Class(sys, alpha, G);
fon_no_input.fsim(u_zero, w_zero, x0_scaled, e_zero, J);

% Define what to estimate
A_known = zeros(n, n);
A2est = ~(A == 0);  % Estimate where A is non-zero
alpha_known = 0.5 * ones(n, 1);
alpha2est = ones(n, 1);  % Estimate all alphas

% Estimate A and alpha
tic;
[est_alpha, est_A, ~] = no_input_estimation_with_known_params(...
    fon_no_input.x, J, A_known, A2est, alpha_known, alpha2est);
t_est = toc;

fprintf('  - Estimation completed in %.2f seconds\n', t_est);
fprintf('  - True alpha: [%.3f, %.3f, %.3f, ...]\n', alpha(1), alpha(2), alpha(3));
fprintf('  - Est. alpha: [%.3f, %.3f, %.3f, ...]\n', est_alpha(1), est_alpha(2), est_alpha(3));

%% Step 5: Parameter Estimation - Stage 2 (B and G)
fprintf('\nStep 5: Estimating B and G matrices...\n');

% Simulate with inputs
fon_input = FON_Graph_Class(sys, alpha, G);
fon_input.fsim(u, w, x0_scaled, e, J);

% Compute fractional derivatives
Z = FON_z_sim(fon_input.x, est_alpha, J);
X = fon_input.x(:, 1:end-1);
R = Z - est_A * X;

% Estimate B and G using least squares
B_mask = (B ~= 0);
G_mask = (G ~= 0);
est_B = zeros(size(B));
est_G = zeros(size(G));

for i = 1:n
    row_inputs = [];
    col_mapping = [];
    
    % Collect B inputs
    for j = 1:m
        if B_mask(i, j)
            row_inputs = [row_inputs; u(j, 1:size(R, 2))];
            col_mapping = [col_mapping; [1, j]];
        end
    end
    
    % Collect G inputs
    for j = 1:v
        if G_mask(i, j)
            row_inputs = [row_inputs; w(j, 1:size(R, 2))];
            col_mapping = [col_mapping; [2, j]];
        end
    end
    
    % Solve if there are inputs
    if ~isempty(row_inputs)
        row_params = R(i, :) / row_inputs;
        for j = 1:size(row_params, 2)
            if col_mapping(j, 1) == 1
                est_B(i, col_mapping(j, 2)) = row_params(j);
            else
                est_G(i, col_mapping(j, 2)) = row_params(j);
            end
        end
    end
end

fprintf('  - B and G matrices estimated\n');

%% Step 6: Validation with New Test Inputs
fprintf('\nStep 6: Validating with new test inputs...\n');

% Generate new test inputs (different from training)
N_test = 30;
u_test = [sin((1:N_test)*2*pi/20); cos((1:N_test)*2*pi/15)];
w_test = 0.2 * randn(v, N_test);
e_test = 0.05 * randn(o, N_test);

% Simulate true system
fon_test_true = FON_Graph_Class(sys, alpha, G);
fon_test_true.fsim(u_test, w_test, x0_scaled, e_test, J);
y_true = fon_test_true.y;

% Simulate with estimated parameters
sys_est = ss(est_A, est_B, C, D);
fon_test_est = FON_Graph_Class(sys_est, est_alpha, est_G);
fon_test_est.fsim(u_test, w_test, x0_scaled, e_test, J);
y_est = fon_test_est.y;

%% Step 7: Compute Performance Metrics
fprintf('\nStep 7: Computing performance metrics...\n');

% Mean Squared Error
mse = mean((y_true - y_est).^2, 2);
mse_avg = mean(mse);

% Correlation
rho = zeros(o, 1);
for i = 1:o
    rho(i) = corr(y_true(i, :)', y_est(i, :)');
end

fprintf('\n========== RESULTS ==========\n');
fprintf('Mean Squared Error per output:\n');
for i = 1:o
    fprintf('  Output %d: %.6f\n', i, mse(i));
end
fprintf('Average MSE: %.6f\n\n', mse_avg);

fprintf('Pearson Correlation per output:\n');
for i = 1:o
    fprintf('  Output %d: %.4f\n', i, rho(i));
end

%% Step 8: Visualize Results
fprintf('\nStep 8: Generating plots...\n');

figure('Position', [100, 100, 1000, 700]);

% Plot each output
for i = 1:o
    subplot(o + 2, 1, i);
    plot(y_true(i, :), 'b-', 'LineWidth', 1.5); hold on;
    plot(y_est(i, :), 'r--', 'LineWidth', 1.5);
    grid on;
    title(sprintf('Output %d (MSE: %.4f, ρ: %.3f)', i, mse(i), rho(i)));
    ylabel(sprintf('y_{%d}[k]', i));
    legend('True', 'Estimated', 'Location', 'best');
end

% Plot controllable inputs
subplot(o + 2, 1, o + 1);
plot(u_test', 'LineWidth', 1.5);
grid on;
title('Test Inputs (u)');
ylabel('u[k]');
legend(arrayfun(@(i) sprintf('u_%d', i), 1:m, 'UniformOutput', false));

% Plot uncontrollable input
subplot(o + 2, 1, o + 2);
plot(w_test', 'LineWidth', 1.5);
grid on;
title('Disturbance (w)');
xlabel('Time Step k');
ylabel('w[k]');

fprintf('\n========================================\n');
fprintf('  Example completed successfully!\n');
fprintf('========================================\n');
