%% ========================================================================
%  FRACTIONAL-ORDER NETWORK (FON) PARAMETER ESTIMATION WITH DECOMPOSITION
%  ========================================================================
%  This code implements the methodology from:
%  "Structural Identifiability in Fractional-Order Networks"
%  
%  Main Steps:
%  1. Generate/load a fractional-order network system
%  2. Decompose the network into subsystems using graph theory
%  3. Estimate parameters for the full system (baseline)
%  4. Estimate parameters for each subsystem independently
%  5. Reconstruct the full system from subsystem estimates
%  6. Validate reconstruction accuracy with new test inputs
%  ========================================================================

clearvars -except;
close all;
clc;

%% ========================================================================
%  SECTION 1: NETWORK GENERATION
%  ========================================================================
%  Generate sparse fractional-order network systems for testing.
%  Multiple trials allow statistical validation of the approach.

% Network dimensions
n = 10;  % Number of state variables
m = 3;   % Number of controllable inputs (u)
v = 2;   % Number of uncontrollable inputs (w) - typically disturbances
o = 3;   % Number of outputs (y)

% Generate multiple random network instances for robust testing
n_trial = 50;  % Number of random networks to test
[A_trial, B_trial, C_trial, D_trial, G_trial, alpha_trial, x0_trial] = ...
    generateSparseNetwork(n, m, v, o, 'high', n_trial);

% Loop through each trial network
for trial_test = 1:n_trial
    
    %% ====================================================================
    %  SECTION 1.1: EXTRACT SYSTEM MATRICES FOR CURRENT TRIAL
    %  ====================================================================
    %  Extract the system matrices that define the fractional-order network:
    %  Δ^α x[k+1] = A*x[k] + B*u[k] + G*w[k]  (with fractional derivative)
    %  y[k] = C*x[k] + D*u[k]
    
    A = A_trial{trial_test, 1};      % State transition matrix (n × n)
    B = B_trial{trial_test, 1};      % Controllable input matrix (n × m)
    C = C_trial{trial_test, 1};      % Output matrix (o × n)
    D = D_trial{trial_test, 1};      % Direct feedthrough matrix (o × m)
    x0 = x0_trial{trial_test, 1};    % Initial state vector (n × 1)
    alpha = alpha_trial{trial_test, 1};  % Fractional orders (n × 1)
    G = G_trial{trial_test, 1};      % Uncontrollable input matrix (n × v)
    
    % Scale initial conditions for better numerical behavior
    x0 = x0 * 100;
    
    % Create state-space system and FON object with graph structure
    sys = ss(A, B, C, D);
    fon = FON_Graph_Class(sys, alpha, G);
    
    %% ====================================================================
    %  SECTION 2: HIERARCHICAL DECOMPOSITION
    %  ====================================================================
    %  Decompose the network into smaller subsystems based on input-output
    %  reachability. This enables identification of subsystems independently,
    %  significantly reducing computational complexity.
    %  
    %  Reference: Algorithm 1 in the paper
    
    sub_systems = fon.hierarchicalDecomposition();
    fprintf('Network decomposed into %d subsystems\n', length(sub_systems));
    
    %% ====================================================================
    %  SECTION 3: FULL SYSTEM PARAMETER ESTIMATION (BASELINE)
    %  ====================================================================
    %  Estimate parameters for the original full network as a baseline.
    %  This uses a two-stage approach:
    %    Stage 1: Estimate A and alpha from zero-input data
    %    Stage 2: Estimate B and G from input-driven data
    
    % Simulation parameters
    N = 100;  % Number of time steps
    J = N;    % Memory length for fractional-order simulation
    
    % ---------------------- STAGE 1: Zero-Input Simulation ----------------
    % Simulate with no controllable input to identify A and alpha
    u_zero = zeros(fon.m, N);
    
    % Generate process and measurement noise
    rng(42);  % Set seed for reproducibility
    noise_sigma_w = 0.5;   % Standard deviation of process noise
    noise_sigma_e = 0.25;  % Standard deviation of measurement noise
    w_noise = noise_sigma_w * randn(fon.v, N);
    e_noise = noise_sigma_e * randn(fon.o, N);
    
    % Simulate fractional-order network with zero controllable input
    fon_no_input = FON_Graph_Class(sys, alpha, G);
    fon_no_input.fsim(u_zero, w_noise, x0, e_noise, J);
    
    % Define structural constraints for A matrix estimation
    % A_known: entries we know exactly (typically zeros from structure)
    % A2est: binary mask (1 = estimate this entry, 0 = use known value)
    A_known = zeros(n, n);   % Known structural zeros
    A2est = ~(A == 0);       % Estimate where A is non-zero
    
    % Define alpha constraints (estimate all for full system)
    alpha_known = 0.5 * ones(n, 1);  % Initial guess
    alpha2est = ones(n, 1);          % Estimate all alpha values
    
    % Estimate A and alpha using zero-input data
    fprintf('\nEstimating full system A and alpha...\n');
    t_start_est = tic;
    [est_alpha, est_A, est_error] = ...
        no_input_estimation_with_known_params(fon_no_input.x, J, A_known, A2est, alpha_known, alpha2est);
    t_end_est = toc(t_start_est);
    fprintf('Full system A estimation completed in %.3f seconds\n', t_end_est);
    
    % ---------------------- STAGE 2: Input-Driven Simulation -------------
    % Simulate with step inputs to identify B and G matrices
    u_step = ones(fon.m, N);  % Step input signal
    
    % Generate new noise realizations for this simulation
    rng(42);  % Same seed for consistency
    w_noise_input = noise_sigma_w * randn(fon.v, N);
    e_noise_input = noise_sigma_e * randn(fon.o, N);
    
    % Simulate with step input
    fon_input = FON_Graph_Class(sys, alpha, G);
    fon_input.fsim(u_step, w_noise_input, x0, e_noise_input, J);
    
    % Calculate Z = fractional derivative of states using estimated alpha
    % Z represents the left-hand side of: Δ^α x[k+1] = A*x[k] + B*u[k] + G*w[k]
    Z = FON_z_sim(fon_input.x, est_alpha, J);
    
    % Calculate residual: R = Z - A*X
    % This residual should equal B*u + G*w
    X = fon_input.x(:, 1:end-1);
    R = Z - est_A * X;
    
    % Create structure masks from the network graph
    % These masks indicate which entries can be non-zero
    B_mask = (B ~= 0);  % 1 where B entries exist in true system
    G_mask = (G ~= 0);  % 1 where G entries exist in true system
    
    % Initialize estimated matrices
    est_B = zeros(size(B));
    est_G = zeros(size(G));
    
    % Estimate B and G row by row using least squares
    % For each state, solve: R(i,:) = B(i,:)*u + G(i,:)*w
    fprintf('Estimating full system B and G...\n');
    for i = 1:size(A, 1)
        % Collect only the relevant inputs for this state
        row_inputs = [];
        col_mapping = [];  % Track which matrix and column each input belongs to
        
        % Add controllable inputs that affect this state
        for j = 1:size(B, 2)
            if B_mask(i, j)
                row_inputs = [row_inputs; u_step(j, 1:size(R, 2))];
                col_mapping = [col_mapping; [1, j]];  % Type 1 = B matrix
            end
        end
        
        % Add uncontrollable inputs that affect this state
        for j = 1:size(G, 2)
            if G_mask(i, j)
                row_inputs = [row_inputs; w_noise_input(j, 1:size(R, 2))];
                col_mapping = [col_mapping; [2, j]];  % Type 2 = G matrix
            end
        end
        
        % Solve least squares if there are inputs to estimate
        if ~isempty(row_inputs)
            row_params = R(i, :) / row_inputs;  % Least-squares solution
            
            % Map results back to B and G matrices
            for j = 1:size(row_params, 2)
                if col_mapping(j, 1) == 1  % B matrix
                    est_B(i, col_mapping(j, 2)) = row_params(j);
                else  % G matrix
                    est_G(i, col_mapping(j, 2)) = row_params(j);
                end
            end
        end
    end
    
    %% ====================================================================
    %  SECTION 4: SUBSYSTEM PARAMETER ESTIMATION
    %  ====================================================================
    %  Estimate parameters for each subsystem independently.
    %  Key innovation: Reuse previously estimated parameters to avoid
    %  redundant computation when subsystems share states.
    
    % Storage for subsystem estimates
    sub_est_alpha = cell(length(sub_systems), 1);
    sub_est_A = cell(length(sub_systems), 1);
    sub_est_B = cell(length(sub_systems), 1);
    sub_est_G = cell(length(sub_systems), 1);
    
    % Create mapping from subsystem states to original system states
    % This allows us to extract the correct data for each subsystem
    state_mapping = cell(length(sub_systems), 1);
    for i = 1:length(sub_systems)
        subsys = sub_systems{i};
        state_mapping{i} = zeros(subsys.n, 1);
        
        % Find corresponding indices in the original system
        for j = 1:subsys.n
            state_id = subsys.G.X(j);  % State ID in subsystem graph
            
            % Locate this state in the original system
            for k = 1:fon.n
                if fon.G.X(k) == state_id
                    state_mapping{i}(j) = k;
                    break;
                end
            end
        end
    end
    
    % Initialize global estimate matrices
    % These accumulate estimates as we process subsystems
    global_A_est = zeros(n, n);
    global_B_est = zeros(n, m);
    global_G_est = zeros(n, v);
    global_alpha_est = zeros(n, 1);  % Track estimated alphas
    
    % Process each subsystem sequentially
    for i = 1:length(sub_systems)
        subsys = sub_systems{i};
        
        % Skip empty subsystems
        if subsys.n == 0
            sub_est_alpha{i} = [];
            sub_est_A{i} = [];
            sub_est_B{i} = zeros(0, max(1, subsys.m));
            sub_est_G{i} = zeros(0, max(1, subsys.v));
            continue;
        end
        
        fprintf('\n========== SUBSYSTEM %d ==========\n', i);
        fprintf('States: %d, Inputs: %d, Disturbances: %d\n', ...
            subsys.n, subsys.m, subsys.v);
        
        % ---------------------- Extract Subsystem State Trajectories -------
        % Pull out the state trajectories that belong to this subsystem
        subsys_states_no_input = zeros(subsys.n, size(fon_no_input.x, 2));
        subsys_states_input = zeros(subsys.n, size(fon_input.x, 2));
        
        for j = 1:subsys.n
            orig_idx = state_mapping{i}(j);
            if orig_idx > 0
                subsys_states_no_input(j, :) = fon_no_input.x(orig_idx, :);
                subsys_states_input(j, :) = fon_input.x(orig_idx, :);
            else
                warning('State mapping failed for subsystem %d, state %d', i, j);
            end
        end
        
        % ---------------------- Build Edge Existence Matrix ---------------
        % Determine which entries in the subsystem's A matrix should be non-zero
        % based on the graph structure
        edge_exists = zeros(subsys.n, subsys.n);
        
        for j = 1:size(subsys.G.E, 1)
            from_node = subsys.G.E(j, 1);
            to_node = subsys.G.E(j, 2);
            
            % Check if this is a state-to-state edge
            if ismember(from_node, subsys.G.X) && ismember(to_node, subsys.G.X)
                from_idx = find(subsys.G.X == from_node);
                to_idx = find(subsys.G.X == to_node);
                
                if ~isempty(from_idx) && ~isempty(to_idx)
                    edge_exists(to_idx, from_idx) = 1;
                end
            end
        end
        
        % ---------------------- Prepare A and Alpha Estimation ------------
        % Reuse previously estimated parameters where possible
        [subsys_A_known, subsys_A2est] = ...
            prepare_subsystem_estimation(state_mapping{i}, global_A_est, edge_exists);
        
        % Prepare alpha masks - reuse from full system or previous subsystems
        subsys_alpha_known = zeros(subsys.n, 1);
        subsys_alpha2est = ones(subsys.n, 1);  % Initially estimate all
        
        for j = 1:subsys.n
            orig_idx = state_mapping{i}(j);
            % Check if we already have an estimate for this alpha
            if global_alpha_est(orig_idx) ~= 0
                % Reuse from previous subsystem estimate
                subsys_alpha_known(j) = global_alpha_est(orig_idx);
                subsys_alpha2est(j) = 0;  % Don't re-estimate
            elseif est_alpha(orig_idx) ~= 0
                % Reuse from full system estimation
                subsys_alpha_known(j) = est_alpha(orig_idx);
                subsys_alpha2est(j) = 0;  % Don't re-estimate
            else
                subsys_alpha_known(j) = 0.5;  % Default starting value
            end
        end
        
        % Report parameter reuse statistics
        n_reused_A = sum(subsys_A2est(:) == 0 & subsys_A_known(:) ~= 0);
        n_to_estimate_A = sum(subsys_A2est(:));
        n_reused_alpha = sum(subsys_alpha2est == 0);
        n_to_estimate_alpha = sum(subsys_alpha2est);
        
        fprintf('A matrix: Reusing %d parameters, estimating %d new parameters\n', ...
            n_reused_A, n_to_estimate_A);
        fprintf('Alpha vector: Reusing %d parameters, estimating %d new parameters\n', ...
            n_reused_alpha, n_to_estimate_alpha);
        
        % Estimate subsystem A and alpha from zero-input data
        t_start_sub = tic;
        [sub_est_alpha{i}, sub_est_A{i}, sub_est_error] = ...
            no_input_estimation_with_known_params(subsys_states_no_input, J, ...
            subsys_A_known, subsys_A2est, subsys_alpha_known, subsys_alpha2est);
        t_end_sub = toc(t_start_sub);
        
        % Track subsystem estimation time for total reconstruction time
        if i == 1
            total_subsystem_time = t_end_sub;  % Initialize on first subsystem
        else
            total_subsystem_time = total_subsystem_time + t_end_sub;
        end
        
        fprintf('Subsystem A estimation completed in %.3f seconds\n', t_end_sub);
        
        % Update global A and alpha estimates with newly estimated parameters
        for row = 1:subsys.n
            for col = 1:subsys.n
                orig_row = state_mapping{i}(row);
                orig_col = state_mapping{i}(col);
                global_A_est(orig_row, orig_col) = sub_est_A{i}(row, col);
            end
        end
        
        for j = 1:subsys.n
            orig_idx = state_mapping{i}(j);
            global_alpha_est(orig_idx) = sub_est_alpha{i}(j);
        end
        
        % ---------------------- Initialize B and G Matrices ---------------
        % Handle edge cases where subsystem has no inputs
        if subsys.m == 0
            sub_est_B{i} = zeros(subsys.n, 1);
        else
            sub_est_B{i} = zeros(subsys.n, subsys.m);
        end
        
        if subsys.v == 0
            sub_est_G{i} = zeros(subsys.n, 1);
        else
            sub_est_G{i} = zeros(subsys.n, subsys.v);
        end
        
        % Skip input estimation if subsystem has no inputs
        if subsys.m == 0 && subsys.v == 0
            continue;
        end
        
        % ---------------------- Map Inputs to Subsystem -------------------
        % Create mappings between subsystem inputs and original system inputs
        u_subsys = zeros(subsys.m, size(u_step, 2));
        w_subsys = zeros(subsys.v, size(w_noise_input, 2));
        
        input_mapping_B = zeros(subsys.m, 1);
        input_mapping_G = zeros(subsys.v, 1);
        
        % Map controllable inputs
        if subsys.m > 0
            for j = 1:subsys.m
                u_id = subsys.G.U(j) - 100;  % Input ID in graph
                for k = 1:fon.m
                    orig_u_id = fon.G.U(k) - 100;
                    if orig_u_id == u_id
                        u_subsys(j, :) = u_step(k, :);
                        input_mapping_B(j) = k;
                        break;
                    end
                end
            end
        end
        
        % Map uncontrollable inputs (disturbances)
        if subsys.v > 0
            for j = 1:subsys.v
                w_id = subsys.G.W(j);
                for k = 1:fon.v
                    if fon.G.W(k) == w_id
                        w_subsys(j, :) = w_noise_input(k, :);
                        input_mapping_G(j) = k;
                        break;
                    end
                end
            end
        end
        
        % ---------------------- Calculate Residual ------------------------
        % Compute Z and residual for this subsystem
        Z_subsys = FON_z_sim(subsys_states_input, sub_est_alpha{i}, J);
        X_subsys = subsys_states_input(:, 1:end-1);
        R_subsys = Z_subsys - sub_est_A{i} * X_subsys;
        
        % ---------------------- Create Structure Masks --------------------
        % Determine which B and G entries should be estimated
        subsys_B_mask = zeros(subsys.n, subsys.m);
        subsys_G_mask = zeros(subsys.n, subsys.v);
        
        % Populate B mask from graph edges
        for j = 1:size(subsys.G.E, 1)
            from_node = subsys.G.E(j, 1);
            to_node = subsys.G.E(j, 2);
            
            % Check for controllable input-to-state edges
            if ismember(from_node, subsys.G.U) && ismember(to_node, subsys.G.X)
                input_idx = find(subsys.G.U == from_node);
                state_idx = find(subsys.G.X == to_node);
                
                if ~isempty(input_idx) && ~isempty(state_idx)
                    subsys_B_mask(state_idx, input_idx) = 1;
                end
            end
        end
        
        % Populate G mask from graph edges
        for j = 1:size(subsys.G.E, 1)
            from_node = subsys.G.E(j, 1);
            to_node = subsys.G.E(j, 2);
            
            % Check for uncontrollable input-to-state edges
            if ismember(from_node, subsys.G.W) && ismember(to_node, subsys.G.X)
                input_idx = find(subsys.G.W == from_node);
                state_idx = find(subsys.G.X == to_node);
                
                if ~isempty(input_idx) && ~isempty(state_idx)
                    subsys_G_mask(state_idx, input_idx) = 1;
                end
            end
        end
        
        % ---------------------- Prepare B and G Estimation ----------------
        % Reuse previously estimated parameters
        [subsys_B_known, subsys_B2est, subsys_G_known, subsys_G2est] = ...
            prepare_subsystem_BG_estimation(state_mapping{i}, input_mapping_B, ...
            input_mapping_G, global_B_est, global_G_est, subsys_B_mask, subsys_G_mask);
        
        % Report parameter reuse statistics
        n_reused_B = sum(subsys_B2est(:) == 0 & subsys_B_known(:) ~= 0);
        n_to_estimate_B = sum(subsys_B2est(:));
        n_reused_G = sum(subsys_G2est(:) == 0 & subsys_G_known(:) ~= 0);
        n_to_estimate_G = sum(subsys_G2est(:));
        fprintf('B matrix: Reusing %d parameters, estimating %d new parameters\n', ...
            n_reused_B, n_to_estimate_B);
        fprintf('G matrix: Reusing %d parameters, estimating %d new parameters\n', ...
            n_reused_G, n_to_estimate_G);
        
        % Copy known values to output matrices
        sub_est_B{i} = subsys_B_known;
        sub_est_G{i} = subsys_G_known;
        
        % ---------------------- Estimate Unknown B and G Entries ----------
        % For each state, estimate only the unknown parameters
        for j = 1:subsys.n
            % Collect inputs that need estimation for this state
            row_inputs = [];
            col_mapping = [];
            
            % Add controllable inputs that need estimation
            for k = 1:subsys.m
                if subsys_B2est(j, k) == 1
                    row_inputs = [row_inputs; u_subsys(k, 1:size(R_subsys, 2))];
                    col_mapping = [col_mapping; [1, k]];
                end
            end
            
            % Add uncontrollable inputs that need estimation
            for k = 1:subsys.v
                if subsys_G2est(j, k) == 1
                    row_inputs = [row_inputs; w_subsys(k, 1:size(R_subsys, 2))];
                    col_mapping = [col_mapping; [2, k]];
                end
            end
            
            % Solve least squares if there are unknowns
            if ~isempty(row_inputs)
                % Adjust residual by subtracting known parameter contributions
                R_adjusted = R_subsys(j, :);
                
                % Subtract known B contributions
                for k = 1:subsys.m
                    if subsys_B2est(j, k) == 0 && subsys_B_known(j, k) ~= 0
                        R_adjusted = R_adjusted - ...
                            subsys_B_known(j, k) * u_subsys(k, 1:size(R_subsys, 2));
                    end
                end
                
                % Subtract known G contributions
                for k = 1:subsys.v
                    if subsys_G2est(j, k) == 0 && subsys_G_known(j, k) ~= 0
                        R_adjusted = R_adjusted - ...
                            subsys_G_known(j, k) * w_subsys(k, 1:size(R_subsys, 2));
                    end
                end
                
                % Solve for unknown parameters
                row_params = R_adjusted / row_inputs;
                
                % Map results back to B and G matrices
                for k = 1:size(row_params, 2)
                    if col_mapping(k, 1) == 1  % B matrix
                        sub_est_B{i}(j, col_mapping(k, 2)) = row_params(k);
                    else  % G matrix
                        sub_est_G{i}(j, col_mapping(k, 2)) = row_params(k);
                    end
                end
            end
        end
        
        % ---------------------- Update Global Estimates -------------------
        % Store estimated B and G parameters in global matrices
        for row = 1:subsys.n
            for col = 1:subsys.m
                if subsys_B_mask(row, col) == 1
                    orig_row = state_mapping{i}(row);
                    orig_col = input_mapping_B(col);
                    global_B_est(orig_row, orig_col) = sub_est_B{i}(row, col);
                end
            end
        end
        
        for row = 1:subsys.n
            for col = 1:subsys.v
                if subsys_G_mask(row, col) == 1
                    orig_row = state_mapping{i}(row);
                    orig_col = input_mapping_G(col);
                    global_G_est(orig_row, orig_col) = sub_est_G{i}(row, col);
                end
            end
        end
    end
    
    %% ====================================================================
    %  SECTION 5: SYSTEM RECONSTRUCTION FROM SUBSYSTEMS
    %  ====================================================================
    %  Reconstruct the full system by combining subsystem estimates.
    %  This process handles duplicate states that may appear in multiple
    %  subsystems and estimates connections between subsystems.
    
    fprintf('\n========== SYSTEM RECONSTRUCTION ==========\n');
    
    % ---------------------- Build Initial Combined System ----------------
    % First, build matrices that include all subsystem states (with duplicates)
    alpha_bar = vertcat(sub_est_alpha{:});
    total_states = length(alpha_bar);
    
    % Initialize B_bar and G_bar
    B_bar = zeros(total_states, fon.m);
    G_bar = zeros(total_states, fon.v);
    
    % Fill B_bar and G_bar by mapping subsystem estimates to full system
    current_row = 0;
    for i = 1:length(sub_systems)
        subsys = sub_systems{i};
        n_states = subsys.n;
        
        if n_states == 0
            continue;
        end
        
        % Map each controllable input
        for input_idx = 1:fon.m
            % Check if this input affects any state in this subsystem
            affects_subsystem = false;
            for j = 1:n_states
                orig_state_idx = state_mapping{i}(j);
                if B(orig_state_idx, input_idx) ~= 0
                    affects_subsystem = true;
                    break;
                end
            end
            
            % If input affects subsystem, find and copy its column
            if affects_subsystem
                subsys_input_col = 0;
                for j = 1:length(subsys.G.U)
                    if subsys.G.U(j) == fon.G.U(input_idx)
                        subsys_input_col = j;
                        break;
                    end
                end
                
                if subsys_input_col > 0
                    B_bar(current_row+1:current_row+n_states, input_idx) = ...
                        sub_est_B{i}(:, subsys_input_col);
                end
            end
        end
        
        % Map each uncontrollable input
        for input_idx = 1:fon.v
            affects_subsystem = false;
            for j = 1:n_states
                orig_state_idx = state_mapping{i}(j);
                if G(orig_state_idx, input_idx) ~= 0
                    affects_subsystem = true;
                    break;
                end
            end
            
            if affects_subsystem
                subsys_w_col = 0;
                for j = 1:length(subsys.G.W)
                    if subsys.G.W(j) == fon.G.W(input_idx)
                        subsys_w_col = j;
                        break;
                    end
                end
                
                if subsys_w_col > 0
                    G_bar(current_row+1:current_row+n_states, input_idx) = ...
                        sub_est_G{i}(:, subsys_w_col);
                end
            end
        end
        
        current_row = current_row + n_states;
    end
    
    % ---------------------- Build Output Matrix C_bar --------------------
    C_subsys = cell(length(sub_systems), 1);
    for i = 1:length(sub_systems)
        subsys = sub_systems{i};
        if subsys.n == 0
            C_subsys{i} = [];
            continue;
        end
        
        C_subsys{i} = zeros(subsys.o, subsys.n);
        
        % Map each output in the subsystem
        for j = 1:subsys.o
            output_id = subsys.G.Y(j);
            
            % Find corresponding output in original system
            for k = 1:fon.o
                if output_id == fon.G.Y(k)
                    % Find which state this output measures
                    for m = 1:size(subsys.G.E, 1)
                        edge = subsys.G.E(m, :);
                        if edge(2) == output_id && ismember(edge(1), subsys.G.X)
                            state_idx = find(subsys.G.X == edge(1));
                            C_subsys{i}(j, state_idx) = 1;
                        end
                    end
                    break;
                end
            end
        end
    end
    
    C_bar = blkdiag(C_subsys{:});
    D_bar = zeros(size(C_bar, 1), size(B_bar, 2));
    
    % ---------------------- Remove Duplicate States ----------------------
    % Identify which states appear multiple times in the combined system
    all_states = [];
    for i = 1:length(sub_systems)
        all_states = [all_states; state_mapping{i}];
    end
    
    fprintf('States before duplicate removal: %d\n', length(all_states));
    
    % Keep only first occurrence of each state
    [unique_states, first_idx] = unique(all_states, 'stable');
    keep_indices = sort(first_idx);
    
    fprintf('States after duplicate removal: %d\n', length(unique_states));
    
    % Remove duplicate rows/columns
    B_bar = B_bar(keep_indices, :);
    G_bar = G_bar(keep_indices, :);
    C_bar = C_bar(:, keep_indices);
    alpha_bar = alpha_bar(keep_indices);
    D_bar = zeros(size(C_bar, 1), size(B_bar, 2));
    
    % ---------------------- Iterative A Matrix Construction --------------
    % Build A_bar iteratively, estimating connections between subsystems
    fprintf('\n--- Building A_bar Iteratively ---\n');
    
    % Start with first subsystem
    A_bar = sub_est_A{1};
    alpha_bar = sub_est_alpha{1};
    x0_bar = x0(state_mapping{1});
    current_size = size(A_bar, 1);
    
    incorporated_states = state_mapping{1};
    remaining_subsystems = 2:length(sub_systems);
    
    fprintf('Starting with subsystem 1: %d states\n', current_size);
    
    % Initialize timing trackers
    total_offdiag_time = 0;
    n_offdiag_estimations = 0;
    offdiag_times = [];  % Store individual times for statistics
    
    % Iteratively add remaining subsystems
    for i = remaining_subsystems
        if isempty(sub_est_A{i}) || isempty(state_mapping{i})
            continue;
        end
        
        fprintf('\n--- Processing Subsystem %d ---\n', i);
        
        % Remove duplicate states from this subsystem
        shared_states = intersect(state_mapping{i}, incorporated_states);
        
        if ~isempty(shared_states)
            fprintf('Shared states: [%s]\n', num2str(shared_states'));
            
            % Find positions of shared states
            shared_positions = [];
            for shared_state = shared_states'
                pos = find(state_mapping{i} == shared_state);
                shared_positions = [shared_positions; pos];
            end
            
            % Keep only non-shared positions
            all_positions = 1:length(state_mapping{i});
            keep_positions = setdiff(all_positions, shared_positions);
            
            fprintf('Removing duplicates at positions [%s], keeping [%s]\n', ...
                num2str(shared_positions'), num2str(keep_positions));
            
            % Reduce subsystem matrices
            A_next = sub_est_A{i}(keep_positions, keep_positions);
            B_next = sub_est_B{i}(keep_positions, :);
            G_next = sub_est_G{i}(keep_positions, :);
            alpha_next = sub_est_alpha{i}(keep_positions);
            state_mapping_next = state_mapping{i}(keep_positions);
            x0_next = x0(state_mapping_next);
        else
            fprintf('No shared states, using full subsystem\n');
            A_next = sub_est_A{i};
            B_next = sub_est_B{i};
            G_next = sub_est_G{i};
            alpha_next = sub_est_alpha{i};
            state_mapping_next = state_mapping{i};
            x0_next = x0(state_mapping_next);
        end
        
        next_size = size(A_next, 1);
        
        if next_size == 0
            fprintf('No new states to add, skipping.\n');
            continue;
        end
        
        fprintf('Adding %d new states\n', next_size);
        
        % Create combined state matrix
        combined_size = current_size + next_size;
        
        % Initialize with known diagonal blocks
        A_known = zeros(combined_size);
        A_known(1:current_size, 1:current_size) = A_bar;
        A_known(current_size+1:end, current_size+1:end) = A_next;
        
        % Determine which off-diagonal entries to estimate
        A2est = zeros(combined_size);
        has_connections = false;
        
        % Check for connections from incorporated states to new states (A21)
        for row = 1:next_size
            orig_row = state_mapping_next(row);
            for col = 1:current_size
                orig_col = incorporated_states(col);
                if A(orig_row, orig_col) ~= 0
                    A2est(current_size + row, col) = 1;
                    has_connections = true;
                end
            end
        end
        
        % Check for connections from new states to incorporated states (A12)
        for row = 1:current_size
            orig_row = incorporated_states(row);
            for col = 1:next_size
                orig_col = state_mapping_next(col);
                if A(orig_row, orig_col) ~= 0
                    A2est(row, current_size + col) = 1;
                    has_connections = true;
                end
            end
        end
        
        % Estimate off-diagonal blocks if connections exist
        if has_connections
            fprintf('Estimating off-diagonal connections...\n');
            
            % Combine state trajectories
            combined_states = [incorporated_states; state_mapping_next];
            combined_states_no_input = zeros(combined_size, size(fon_no_input.x, 2));
            
            for j = 1:combined_size
                orig_idx = combined_states(j);
                combined_states_no_input(j, :) = fon_no_input.x(orig_idx, :);
            end
            
            % Use known alphas from subsystem estimates
            combined_alpha = [alpha_bar; alpha_next];
            combined_alpha2est = zeros(combined_size, 1);  % All alphas are KNOWN!
            
            % Count parameters
            n_unknown_A = sum(A2est(:));
            n_known_diag = sum(A_known(:) ~= 0);
            fprintf('Matrix structure: %d diagonal parameters (known), %d off-diagonal parameters (estimating)\n', ...
                n_known_diag, n_unknown_A);
            
            % Estimate with alpha fixed - only optimize A off-diagonals
            t_start_offdiag = tic;
            [~, A_bar, ~] = no_input_estimation_with_known_params(...
                combined_states_no_input, J, A_known, A2est, combined_alpha, combined_alpha2est);
            t_offdiag = toc(t_start_offdiag);
            
            % Track timing statistics
            total_offdiag_time = total_offdiag_time + t_offdiag;
            n_offdiag_estimations = n_offdiag_estimations + 1;
            offdiag_times = [offdiag_times; t_offdiag];
            
            fprintf('Off-diagonal block estimation completed in %.3f seconds\n', t_offdiag);
        else
            fprintf('No connections found, using block diagonal structure\n');
            A_bar = A_known;
        end
        
        % Update combined system
        alpha_bar = [alpha_bar; alpha_next];
        x0_bar = [x0_bar; x0_next];
        current_size = combined_size;
        incorporated_states = [incorporated_states; state_mapping_next];
    end
    
    % Build final C_bar based on incorporated states
    C_bar = zeros(fon.o, length(incorporated_states));
    for i = 1:fon.o
        state_measured = find(C(i, :) ~= 0);
        if ~isempty(state_measured)
            pos_in_incorporated = find(incorporated_states == state_measured(1));
            if ~isempty(pos_in_incorporated)
                C_bar(i, pos_in_incorporated) = C(i, state_measured(1));
            end
        end
    end
    
    fprintf('\n=== Reconstruction Complete ===\n');
    fprintf('Final system size: %d states\n', length(alpha_bar));
    
    % Report detailed timing statistics
    if n_offdiag_estimations > 0
        fprintf('\n--- Off-Diagonal Block Estimation Timing ---\n');
        fprintf('Number of off-diagonal estimations: %d\n', n_offdiag_estimations);
        fprintf('Total time for off-diagonal blocks: %.3f seconds\n', total_offdiag_time);
        fprintf('Average time per off-diagonal block: %.3f seconds\n', total_offdiag_time / n_offdiag_estimations);
        fprintf('Min time: %.3f seconds, Max time: %.3f seconds\n', min(offdiag_times), max(offdiag_times));
    else
        fprintf('\n--- Off-Diagonal Block Estimation Timing ---\n');
        fprintf('No off-diagonal blocks required estimation (fully decoupled subsystems)\n');
    end
    
    % Calculate total reconstruction time
    total_reconstruction_time = total_subsystem_time + total_offdiag_time;
    
    fprintf('\n--- Computational Efficiency Comparison ---\n');
    fprintf('Full system A estimation: %.3f seconds\n', t_end_est);
    fprintf('\nSubsystem-based reconstruction breakdown:\n');
    fprintf('  - Subsystem estimations: %.3f seconds\n', total_subsystem_time);
    fprintf('  - Off-diagonal blocks: %.3f seconds\n', total_offdiag_time);
    fprintf('  - Total reconstruction: %.3f seconds\n', total_reconstruction_time);
    
    if total_reconstruction_time < t_end_est
        speedup = t_end_est / total_reconstruction_time;
        fprintf('\nSpeedup factor: %.2fx faster\n', speedup);
        fprintf('Time saved: %.3f seconds (%.1f%% reduction)\n', ...
            t_end_est - total_reconstruction_time, ...
            100 * (t_end_est - total_reconstruction_time) / t_end_est);
    else
        slowdown = total_reconstruction_time / t_end_est;
        fprintf('\nSlowdown factor: %.2fx slower\n', slowdown);
        fprintf('Additional time: %.3f seconds\n', total_reconstruction_time - t_end_est);
    end
    
    %% ====================================================================
    %  SECTION 6: VALIDATION WITH NEW TEST INPUTS
    %  ====================================================================
    %  Validate the reconstruction by comparing outputs from:
    %    1. True original system
    %    2. Directly estimated original system
    %    3. Reconstructed system from subsystem estimates
    
    fprintf('\n========== VALIDATION ==========\n');
    
    % Generate new test input sequences (different from training)
    N_test = 100;
    J_test = N_test;
    
    % Create sinusoidal test inputs
    u_test = zeros(fon.m, N_test);
    u_test(1, :) = 1 + sin(0.2 * (1:N_test)) + noise_sigma_w * randn(1, N_test);
    u_test(2, :) = 2 + cos(0.15 * (1:N_test)) + noise_sigma_w * randn(1, N_test);
    
    % Generate test noise
    rng(100);  % Different seed for validation
    w_test = noise_sigma_w * randn(fon.v, N_test);
    e_test = noise_sigma_e * randn(fon.o, N_test);
    
    % Use same inputs for reconstructed system
    u_bar = u_test;
    w_bar = w_test;
    
    % Setup noise for reconstructed system outputs
    e_bar = zeros(size(C_bar, 1), N_test);
    if size(e_bar, 1) == size(e_test, 1)
        e_bar = e_test;
    else
        % Replicate noise if dimensions differ
        for i = 1:size(e_bar, 1)
            e_bar(i, :) = e_test(min(i, size(e_test, 1)), :);
        end
    end
    
    % ---------------------- Simulate All Three Systems -------------------
    % 1. True system
    fon_true = FON_Graph_Class(sys, alpha, G);
    fon_true.fsim(u_test, w_test, x0, e_test, J_test);
    y_true = fon_true.y;
    
    % 2. Directly estimated system
    sys_est = ss(est_A, est_B, C, D);
    fon_est = FON_Graph_Class(sys_est, est_alpha, est_G);
    fon_est.fsim(u_test, w_test, x0, e_test, J_test);
    y_est = fon_est.y;
    
    % 3. Reconstructed system from subsystems
    sys_bar = ss(A_bar, B_bar, C_bar, D_bar);
    fon_bar = FON_Graph_Class(sys_bar, alpha_bar, G_bar);
    fon_bar.fsim(u_bar, w_bar, x0_bar, e_bar, J_test);
    y_bar = fon_bar.y;
    
    %% ====================================================================
    %  SECTION 7: RESULTS VISUALIZATION AND METRICS
    %  ====================================================================
    
    % Create comprehensive comparison plot
    figure('Position', [100, 100, 1200, 1000]);
    ax = [];
    
    % Color palette (color-blind friendly)
    c_true = [0.0, 0.45, 0.70];   % Blue
    c_est  = [0.90, 0.60, 0.00];  % Orange
    c_bar  = [0.00, 0.62, 0.45];  % Green
    c_u1   = [0.80, 0.47, 0.65];  % Purple
    c_u2   = [0.35, 0.70, 0.90];  % Sky blue
    c_w1   = [0.94, 0.90, 0.25];  % Yellow
    
    % Plot each output
    for i = 1:fon.o
        ax(i) = subplot(fon.o + 2, 1, i);
        plot(y_true(i, :), '-', 'Color', c_true, 'LineWidth', 1.5); hold on;
        plot(y_est(i, :), '--', 'Color', c_est, 'LineWidth', 1.5);
        plot(y_bar(i, :), '-.', 'Color', c_bar, 'LineWidth', 1.5);
        grid on;
        title(sprintf('Output %d', i), 'FontSize', 14, 'FontWeight', 'bold');
        legend('True System', 'Direct Estimation', 'Subsystem Reconstruction', ...
            'FontSize', 10, 'Location', 'best');
        ylabel(sprintf('y_{%d}[k]', i), 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'FontSize', 11);
    end
    
    % Plot controllable inputs
    ax(fon.o + 1) = subplot(fon.o + 2, 1, fon.o + 1);
    plot(u_test(1, :), '-', 'Color', c_u1, 'LineWidth', 1.5); hold on;
    plot(u_test(2, :), '-', 'Color', c_u2, 'LineWidth', 1.5);
    grid on;
    title('Controllable Inputs', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('u[k]', 'FontSize', 12, 'FontWeight', 'bold');
    legend('u_1[k]', 'u_2[k]', 'FontSize', 10, 'Location', 'best');
    set(gca, 'FontSize', 11);
    
    % Plot uncontrollable input
    ax(fon.o + 2) = subplot(fon.o + 2, 1, fon.o + 2);
    plot(w_test(1, :), '-', 'Color', c_w1, 'LineWidth', 1.5);
    grid on;
    title('Uncontrollable Input (Disturbance)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Time Step k', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('w_1[k]', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'FontSize', 11);
    
    % Link all x-axes
    linkaxes(ax, 'x');
    
    % ---------------------- Compute Error Metrics ------------------------
    % Mean Squared Error
    mse_est = mean((y_true - y_est).^2, 2);
    mse_bar = mean((y_true - y_bar).^2, 2);
    
    % Pearson Correlation Coefficient
    rho_est = zeros(fon.o, 1);
    rho_bar = zeros(fon.o, 1);
    for i = 1:fon.o
        rho_est(i) = corr(y_true(i, :)', y_est(i, :)');
        rho_bar(i) = corr(y_true(i, :)', y_bar(i, :)');
    end
    
    % Display results
    fprintf('\n========== PERFORMANCE METRICS ==========\n');
    fprintf('\nMean Squared Error (MSE):\n');
    for i = 1:fon.o
        fprintf('Output %d - Direct Est: %.6f, Subsys Recon: %.6f\n', ...
            i, mse_est(i), mse_bar(i));
    end
    
    fprintf('\nAverage MSE across all outputs:\n');
    fprintf('Direct Estimation: %.6f\n', mean(mse_est));
    fprintf('Subsystem Reconstruction: %.6f\n', mean(mse_bar));
    
    fprintf('\nPearson Correlation Coefficients:\n');
    for i = 1:fon.o
        fprintf('Output %d - Direct Est: %.4f, Subsys Recon: %.4f\n', ...
            i, rho_est(i), rho_bar(i));
    end
    fprintf('\n');
    
end  % End of trial loop

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [subsys_A_known, subsys_A2est] = prepare_subsystem_estimation(...
    state_mapping_i, global_A_est, edge_exists)
% PREPARE_SUBSYSTEM_ESTIMATION - Prepare A matrix estimation for a subsystem
%
% This function determines which entries of the subsystem's A matrix should
% be re-estimated and which can be reused from previous estimates.
%
% Inputs:
%   state_mapping_i - Mapping from subsystem states to original system (n_sub × 1)
%   global_A_est    - Current estimate of full A matrix (n × n)
%   edge_exists     - Binary matrix indicating graph edges (n_sub × n_sub)
%
% Outputs:
%   subsys_A_known  - Known A values for this subsystem (n_sub × n_sub)
%   subsys_A2est    - Binary mask: 1 = estimate, 0 = use known (n_sub × n_sub)

n_subsys = length(state_mapping_i);

% Initialize: start by assuming we need to estimate where edges exist
subsys_A_known = zeros(n_subsys, n_subsys);
subsys_A2est = edge_exists;

% Allow diagonal entries (needed for fractional-order terms)
for i = 1:n_subsys
    subsys_A2est(i, i) = 1;
end

% Check each entry to see if we can reuse a previous estimate
for i = 1:n_subsys
    for j = 1:n_subsys
        % Map to original system indices
        orig_i = state_mapping_i(i);
        orig_j = state_mapping_i(j);
        
        % If this parameter was already estimated, reuse it
        if global_A_est(orig_i, orig_j) ~= 0
            subsys_A_known(i, j) = global_A_est(orig_i, orig_j);
            subsys_A2est(i, j) = 0;  % Don't re-estimate
        end
    end
end
end

function [subsys_B_known, subsys_B2est, subsys_G_known, subsys_G2est] = ...
    prepare_subsystem_BG_estimation(state_mapping_i, input_mapping_B, ...
    input_mapping_G, global_B_est, global_G_est, subsys_B_mask, subsys_G_mask)
% PREPARE_SUBSYSTEM_BG_ESTIMATION - Prepare B and G matrix estimation
%
% This function determines which entries of B and G matrices should be
% re-estimated and which can be reused from previous subsystem estimates.
%
% Inputs:
%   state_mapping_i  - Subsystem to original state mapping (n_sub × 1)
%   input_mapping_B  - Subsystem to original input mapping for B (m_sub × 1)
%   input_mapping_G  - Subsystem to original input mapping for G (v_sub × 1)
%   global_B_est     - Current global B estimate (n × m)
%   global_G_est     - Current global G estimate (n × v)
%   subsys_B_mask    - Binary mask for B structure (n_sub × m_sub)
%   subsys_G_mask    - Binary mask for G structure (n_sub × v_sub)
%
% Outputs:
%   subsys_B_known   - Known B values (n_sub × m_sub)
%   subsys_B2est     - Binary mask for B: 1 = estimate, 0 = known
%   subsys_G_known   - Known G values (n_sub × v_sub)
%   subsys_G2est     - Binary mask for G: 1 = estimate, 0 = known

n_subsys = length(state_mapping_i);
m_subsys = length(input_mapping_B);
v_subsys = length(input_mapping_G);

% Initialize B matrices
subsys_B_known = zeros(n_subsys, m_subsys);
subsys_B2est = subsys_B_mask;  % Start with structure mask

% Initialize G matrices
subsys_G_known = zeros(n_subsys, v_subsys);
subsys_G2est = subsys_G_mask;  % Start with structure mask

% Check each B entry for reusable estimates
for i = 1:n_subsys
    for j = 1:m_subsys
        if subsys_B_mask(i, j) == 1  % Only where connection exists
            % Map to original indices
            orig_i = state_mapping_i(i);
            orig_j = input_mapping_B(j);
            
            % Reuse if already estimated
            if global_B_est(orig_i, orig_j) ~= 0
                subsys_B_known(i, j) = global_B_est(orig_i, orig_j);
                subsys_B2est(i, j) = 0;  % Don't re-estimate
            end
        end
    end
end

% Check each G entry for reusable estimates
for i = 1:n_subsys
    for j = 1:v_subsys
        if subsys_G_mask(i, j) == 1  % Only where connection exists
            % Map to original indices
            orig_i = state_mapping_i(i);
            orig_j = input_mapping_G(j);
            
            % Reuse if already estimated
            if global_G_est(orig_i, orig_j) ~= 0
                subsys_G_known(i, j) = global_G_est(orig_i, orig_j);
                subsys_G2est(i, j) = 0;  % Don't re-estimate
            end
        end
    end
end
end
