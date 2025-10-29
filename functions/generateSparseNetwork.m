% function [A, B, C, D, G, alpha, x0] = generateSparseNetwork(n, m, v, o, sparsity_level)
%     % GENERATESPARSE NETWORK - Creates a sparse network with guaranteed I/O reachability
%     %
%     % Inputs:
%     %   n - number of states
%     %   m - number of controllable inputs
%     %   v - number of uncontrollable inputs
%     %   o - number of outputs
%     %   sparsity_level - (optional) 'low', 'medium', 'high' (default: 'medium')
%     %
%     % Outputs:
%     %   A - State transition matrix (n x n)
%     %   B - Input matrix (n x m) - EACH INPUT AFFECTS EXACTLY ONE STATE
%     %   C - Output matrix (o x n) - EACH OUTPUT MEASURES EXACTLY ONE STATE
%     %   D - Feedthrough matrix (o x m)
%     %   G - Disturbance matrix (n x v)
%     %   alpha - Fractional orders (n x 1)
%     %   x0 - Initial conditions (n x 1)
% 
%     if nargin < 5
%         sparsity_level = 'medium';
%     end
% 
%     % Set sparsity parameters
%     switch sparsity_level
%         case 'low'
%             min_states_per_subnet = 3;
%             max_states_per_subnet = 5;
%             internal_connection_prob = 0.4;
%         case 'medium'
%             min_states_per_subnet = 2;
%             max_states_per_subnet = 4;
%             internal_connection_prob = 0.3;
%         case 'high'
%             min_states_per_subnet = 2;
%             max_states_per_subnet = 3;
%             internal_connection_prob = 0.2;
%         otherwise
%             error('sparsity_level must be low, medium, or high');
%     end
% 
%     % Keep generating until valid
%     all_reachable = false;
%     attempt = 0;
%     max_attempts = 100;
% 
%     while ~all_reachable && attempt < max_attempts
%         attempt = attempt + 1;
% 
%         %% Step 1: Decide how many subnetworks to create
%         % Base number on outputs and states available
%         num_subnetworks = max(ceil(o * 0.8), min(o, ceil(n / max_states_per_subnet)));
%         num_subnetworks = max(2, min(num_subnetworks, n)); % At least 2, at most n
% 
%         %% Step 2: Distribute states to each subnetwork
%         states_per_subnet = distribute_states(n, num_subnetworks, min_states_per_subnet, max_states_per_subnet);
% 
%         % Build subnetwork state indices
%         subnetworks = cell(num_subnetworks, 1);
%         current_state = 1;
%         for i = 1:num_subnetworks
%             num_states = states_per_subnet(i);
%             subnetworks{i}.states = current_state:(current_state + num_states - 1);
%             subnetworks{i}.num_states = num_states;
%             current_state = current_state + num_states;
%         end
% 
%         %% Step 3: Connect states within each subnetwork (toward output)
%         A = zeros(n, n);
% 
%         for i = 1:num_subnetworks
%             states = subnetworks{i}.states;
%             num_states = length(states);
% 
%             if num_states == 1
%                 % Single state - just add self-loop later
%                 subnetworks{i}.initial_node = states(1);
%                 subnetworks{i}.final_node = states(1);
%             else
%                 % Create chain: initial -> ... -> final
%                 subnetworks{i}.initial_node = states(1);
%                 subnetworks{i}.final_node = states(end);
% 
%                 % Main forward chain
%                 for j = 1:(num_states-1)
%                     A(states(j+1), states(j)) = rand() * 0.3 + 0.1;
%                 end
% 
%                 % Add some additional internal connections (sparse)
%                 if num_states > 2
%                     num_internal = randi([0, floor(num_states * internal_connection_prob)]);
%                     for k = 1:num_internal
%                         source = states(randi(num_states-1)); % Not the last state
%                         target = states(randi([2, num_states])); % Not the first state
%                         if source ~= target && A(target, source) == 0
%                             A(target, source) = rand() * 0.3 + 0.1;
%                         end
%                     end
%                 end
%             end
%         end
% 
%         %% Step 4: Connect outputs - one per subnetwork, cycle if more outputs
%         C = zeros(o, n);
%         for i = 1:o
%             subnet_idx = mod(i-1, num_subnetworks) + 1;
%             output_state = subnetworks{subnet_idx}.final_node;
%             C(i, output_state) = 1.0;
%         end
% 
%         %% Step 5: Connect inputs - one per subnetwork (at initial node), cycle if more inputs
%         B = zeros(n, m);
%         for i = 1:m
%             subnet_idx = mod(i-1, num_subnetworks) + 1;
%             input_state = subnetworks{subnet_idx}.initial_node;
% 
%             % Each input affects exactly one state
%             % If this state already has an input, find another state in the same subnet
%             if B(input_state, :) * B(input_state, :)' > 0
%                 % Initial node occupied, try other states in subnet
%                 available_states = subnetworks{subnet_idx}.states;
%                 for candidate_state = available_states
%                     if B(candidate_state, :) * B(candidate_state, :)' == 0
%                         input_state = candidate_state;
%                         break;
%                     end
%                 end
%             end
% 
%             B(input_state, i) = rand() * 0.5 + 0.5;
%         end
% 
%         %% Step 6: Add self-loops for stability
%         for i = 1:n
%             A(i, i) = -(rand() * 0.3 + 0.4);
%         end
% 
%         %% Step 7: Add inter-subnetwork connections (chain structure)
%         % Start from first subnet, create connection to another subnet
%         % Then from that subnet to another, etc., until only one subnet remains
% 
%         available_subnets = 2:num_subnetworks; % Subnets that haven't been connected yet
%         current_subnet = 1; % Start from first subnet
% 
%         while ~isempty(available_subnets)
%             % Pick a random target from available subnets
%             target_idx = randi(length(available_subnets));
%             target_subnet = available_subnets(target_idx);
% 
%             % Create connection from current_subnet to target_subnet
%             % From final node of source to a random node in target
%             source_state = subnetworks{current_subnet}.final_node;
%             target_states = subnetworks{target_subnet}.states;
%             target_state = target_states(randi(length(target_states)));
% 
%             A(target_state, source_state) = rand() * 0.3 + 0.1;
% 
%             % Remove target_subnet from available list
%             available_subnets(target_idx) = [];
% 
%             % Move to target_subnet for next iteration
%             current_subnet = target_subnet;
%         end
%         % Last subnet (current_subnet) has no outgoing inter-subnet connection
% 
%         %% Step 8: Add uncontrollable inputs (disturbances)
%         G = zeros(n, v);
%         for j = 1:v
%             % Each disturbance affects 2-3 states across different subnetworks
%             num_affected = randi([2, min(3, n)]);
%             affected_states = randperm(n, num_affected);
% 
%             for state_idx = affected_states
%                 G(state_idx, j) = rand() * 0.3 + 0.2;
%             end
%         end
% 
%         %% Step 9: Generate fractional orders and initial conditions
%         alpha = rand(n, 1) * 0.6 + 0.2;  % Between 0.2 and 0.8
%         x0 = randn(n, 1) * 0.5;
% 
%         %% Step 10: Verify reachability
%         all_reachable = true;
%         for i = 1:o
%             output_state = find(C(i, :) ~= 0);
%             if isempty(output_state)
%                 all_reachable = false;
%                 break;
%             end
% 
%             reaches_input = checkReachability(A, B, G, output_state(1));
% 
%             if ~reaches_input
%                 all_reachable = false;
%                 break;
%             end
%         end
% 
%         if ~all_reachable && attempt < max_attempts
%             continue;
%         end
%     end
% 
%     % Check if we hit max attempts
%     if attempt >= max_attempts && ~all_reachable
%         error('Failed to generate a valid network after %d attempts. Try adjusting network parameters.', max_attempts);
%     end
% 
%     % Create D matrix
%     D = zeros(o, m);
% end
% 
% function states_distribution = distribute_states(n, num_subnetworks, min_states, max_states)
%     % Distribute n states among num_subnetworks
%     % Each subnetwork gets between min_states and max_states
% 
%     states_distribution = zeros(num_subnetworks, 1);
%     remaining_states = n;
% 
%     for i = 1:num_subnetworks-1
%         % Calculate how many states this subnet can have
%         min_allowed = min_states;
%         max_allowed = min(max_states, remaining_states - (num_subnetworks - i) * min_states);
% 
%         if max_allowed < min_allowed
%             max_allowed = min_allowed;
%         end
% 
%         % Assign random number of states within bounds
%         states_distribution(i) = randi([min_allowed, max_allowed]);
%         remaining_states = remaining_states - states_distribution(i);
%     end
% 
%     % Last subnetwork gets remaining states
%     states_distribution(num_subnetworks) = remaining_states;
% 
%     % Ensure last subnet has at least 1 state
%     if states_distribution(num_subnetworks) < 1
%         % Steal from previous subnetworks
%         for i = (num_subnetworks-1):-1:1
%             if states_distribution(i) > min_states
%                 transfer = min(states_distribution(i) - min_states, 1 - states_distribution(num_subnetworks));
%                 states_distribution(i) = states_distribution(i) - transfer;
%                 states_distribution(num_subnetworks) = states_distribution(num_subnetworks) + transfer;
%                 if states_distribution(num_subnetworks) >= 1
%                     break;
%                 end
%             end
%         end
%     end
% end
% 
% function reaches_input = checkReachability(A, B, G, start_state)
%     % Quick backward reachability check
%     reachable_states = start_state;
%     to_explore = start_state;
%     explored = [];
% 
%     max_iterations = size(A, 1);
%     iteration = 0;
% 
%     while ~isempty(to_explore) && iteration < max_iterations
%         iteration = iteration + 1;
%         current = to_explore(1);
%         to_explore(1) = [];
% 
%         if ismember(current, explored)
%             continue;
%         end
%         explored = [explored, current];
% 
%         % Find predecessors
%         predecessors = find(A(:, current) ~= 0);
% 
%         for pred = predecessors'
%             if ~ismember(pred, reachable_states)
%                 reachable_states = [reachable_states, pred];
%                 to_explore = [to_explore, pred];
%             end
%         end
% 
%         % Check bidirectional
%         successors = find(A(current, :) ~= 0);
%         for succ = successors
%             if A(succ, current) ~= 0 && ~ismember(succ, reachable_states)
%                 reachable_states = [reachable_states, succ];
%                 to_explore = [to_explore, succ];
%             end
%         end
%     end
% 
%     % Check if any input reaches these states
%     reaches_input = false;
%     for state = reachable_states
%         if any(B(state, :) ~= 0) || any(G(state, :) ~= 0)
%             reaches_input = true;
%             break;
%         end
%     end
% end

function [A, B, C, D, G, alpha, x0] = generateSparseNetwork(n, m, v, o, sparsity_level, num_networks)
    % GENERATESPARSE NETWORK - Creates sparse network(s) with guaranteed I/O reachability
    %
    % Inputs:
    %   n - number of states
    %   m - number of controllable inputs
    %   v - number of uncontrollable inputs
    %   o - number of outputs
    %   sparsity_level - (optional) 'low', 'medium', 'high' (default: 'medium')
    %   num_networks - (optional) number of different networks to generate (default: 1)
    %
    % Outputs (if num_networks == 1):
    %   A, B, C, D, G, alpha, x0 - Single matrices/vectors
    %
    % Outputs (if num_networks > 1):
    %   A, B, C, D, G, alpha, x0 - Cell arrays containing num_networks different networks
    
    if nargin < 5
        sparsity_level = 'medium';
    end
    
    if nargin < 6
        num_networks = 1;
    end
    
    % Set sparsity parameters
    switch sparsity_level
        case 'low'
            min_states_per_subnet = 3;
            max_states_per_subnet = 5;
            internal_connection_prob = 0.4;
        case 'medium'
            min_states_per_subnet = 2;
            max_states_per_subnet = 4;
            internal_connection_prob = 0.3;
        case 'high'
            min_states_per_subnet = 2;
            max_states_per_subnet = 3;
            internal_connection_prob = 0.2;
        otherwise
            error('sparsity_level must be low, medium, or high');
    end
    
    % Initialize storage for multiple networks
    if num_networks > 1
        A_all = cell(num_networks, 1);
        B_all = cell(num_networks, 1);
        C_all = cell(num_networks, 1);
        D_all = cell(num_networks, 1);
        G_all = cell(num_networks, 1);
        alpha_all = cell(num_networks, 1);
        x0_all = cell(num_networks, 1);
    end
    
    % Generate num_networks different networks
    for net_idx = 1:num_networks
        
        % Seed randomness differently for each network
        rng('shuffle');  % Use current time for randomness
        pause(0.01);  % Small pause to ensure different seed
        
        % Keep generating until valid and unique
        all_reachable = false;
        attempt = 0;
        max_attempts = 100;
        
        while ~all_reachable && attempt < max_attempts
            attempt = attempt + 1;
            
            %% Step 1: Decide how many subnetworks to create
            % Base number on outputs and states available
            num_subnetworks = max(ceil(o * 0.8), min(o, ceil(n / max_states_per_subnet)));
            num_subnetworks = max(2, min(num_subnetworks, n)); % At least 2, at most n
            
            %% Step 2: Distribute states to each subnetwork
            states_per_subnet = distribute_states(n, num_subnetworks, min_states_per_subnet, max_states_per_subnet);
            
            % Build subnetwork state indices
            subnetworks = cell(num_subnetworks, 1);
            current_state = 1;
            for i = 1:num_subnetworks
                num_states = states_per_subnet(i);
                subnetworks{i}.states = current_state:(current_state + num_states - 1);
                subnetworks{i}.num_states = num_states;
                current_state = current_state + num_states;
            end
            
            %% Step 3: Connect states within each subnetwork (toward output)
            A_temp = zeros(n, n);
            
            for i = 1:num_subnetworks
                states = subnetworks{i}.states;
                num_states = length(states);
                
                if num_states == 1
                    % Single state - just add self-loop later
                    subnetworks{i}.initial_node = states(1);
                    subnetworks{i}.final_node = states(1);
                else
                    % Create chain: initial -> ... -> final
                    subnetworks{i}.initial_node = states(1);
                    subnetworks{i}.final_node = states(end);
                    
                    % Main forward chain
                    for j = 1:(num_states-1)
                        A_temp(states(j+1), states(j)) = rand() * 0.3 + 0.1;
                    end
                    
                    % Add some additional internal connections (sparse)
                    if num_states > 2
                        num_internal = randi([0, floor(num_states * internal_connection_prob)]);
                        for k = 1:num_internal
                            source = states(randi(num_states-1)); % Not the last state
                            target = states(randi([2, num_states])); % Not the first state
                            if source ~= target && A_temp(target, source) == 0
                                A_temp(target, source) = rand() * 0.3 + 0.1;
                            end
                        end
                    end
                end
            end
            
            %% Step 4: Connect outputs - one per subnetwork, cycle if more outputs
            C_temp = zeros(o, n);
            for i = 1:o
                subnet_idx = mod(i-1, num_subnetworks) + 1;
                output_state = subnetworks{subnet_idx}.final_node;
                C_temp(i, output_state) = 1.0;
            end
            
            %% Step 5: Connect inputs - one per subnetwork (at initial node), cycle if more inputs
            B_temp = zeros(n, m);
            for i = 1:m
                subnet_idx = mod(i-1, num_subnetworks) + 1;
                input_state = subnetworks{subnet_idx}.initial_node;
                
                % Each input affects exactly one state
                % If this state already has an input, find another state in the same subnet
                if B_temp(input_state, :) * B_temp(input_state, :)' > 0
                    % Initial node occupied, try other states in subnet
                    available_states = subnetworks{subnet_idx}.states;
                    for candidate_state = available_states
                        if B_temp(candidate_state, :) * B_temp(candidate_state, :)' == 0
                            input_state = candidate_state;
                            break;
                        end
                    end
                end
                
                B_temp(input_state, i) = rand() * 0.5 + 0.5;
            end
            
            %% Step 6: Add self-loops for stability
            for i = 1:n
                A_temp(i, i) = -(rand() * 0.3 + 0.4);
            end
            
            %% Step 7: Add inter-subnetwork connections (chain structure)
            % Start from first subnet, create connection to another subnet
            % Then from that subnet to another, etc., until only one subnet remains
            
            available_subnets = 2:num_subnetworks; % Subnets that haven't been connected yet
            current_subnet = 1; % Start from first subnet
            
            while ~isempty(available_subnets)
                % Pick a random target from available subnets
                target_idx = randi(length(available_subnets));
                target_subnet = available_subnets(target_idx);
                
                % Create connection from current_subnet to target_subnet
                % From final node of source to a random node in target
                source_state = subnetworks{current_subnet}.final_node;
                target_states = subnetworks{target_subnet}.states;
                target_state = target_states(randi(length(target_states)));
                
                A_temp(target_state, source_state) = rand() * 0.3 + 0.1;
                
                % Remove target_subnet from available list
                available_subnets(target_idx) = [];
                
                % Move to target_subnet for next iteration
                current_subnet = target_subnet;
            end
            % Last subnet (current_subnet) has no outgoing inter-subnet connection
            
            %% Step 8: Add uncontrollable inputs (disturbances)
            G_temp = zeros(n, v);
            for j = 1:v
                % Each disturbance affects 2-3 states across different subnetworks
                num_affected = randi([2, min(3, n)]);
                affected_states = randperm(n, num_affected);
                
                for state_idx = affected_states
                    G_temp(state_idx, j) = rand() * 0.3 + 0.2;
                end
            end
            
            %% Step 9: Generate fractional orders and initial conditions
            alpha_temp = rand(n, 1) * 0.6 + 0.2;  % Between 0.2 and 0.8
            x0_temp = randn(n, 1) * 0.5;
            
            %% Step 10: Verify reachability
            all_reachable = true;
            for i = 1:o
                output_state = find(C_temp(i, :) ~= 0);
                if isempty(output_state)
                    all_reachable = false;
                    break;
                end
                
                reaches_input = checkReachability(A_temp, B_temp, G_temp, output_state(1));
                
                if ~reaches_input
                    all_reachable = false;
                    break;
                end
            end
            
            %% Step 11: Check uniqueness (if generating multiple networks)
            if all_reachable && num_networks > 1 && net_idx > 1
                is_unique = true;
                for prev_idx = 1:(net_idx-1)
                    % Check if this network is too similar to previous ones
                    if isSimilarNetwork(A_temp, B_temp, C_temp, G_temp, ...
                                       A_all{prev_idx}, B_all{prev_idx}, C_all{prev_idx}, G_all{prev_idx})
                        is_unique = false;
                        break;
                    end
                end
                if ~is_unique
                    all_reachable = false;  % Force regeneration
                end
            end
            
            if ~all_reachable && attempt < max_attempts
                continue;
            end
        end
        
        % Check if we hit max attempts
        if attempt >= max_attempts && ~all_reachable
            error('Failed to generate a valid network %d after %d attempts. Try adjusting network parameters.', net_idx, max_attempts);
        end
        
        % Create D matrix
        D_temp = zeros(o, m);
        
        % Store the generated network
        if num_networks > 1
            A_all{net_idx} = A_temp;
            B_all{net_idx} = B_temp;
            C_all{net_idx} = C_temp;
            D_all{net_idx} = D_temp;
            G_all{net_idx} = G_temp;
            alpha_all{net_idx} = alpha_temp;
            x0_all{net_idx} = x0_temp;
        else
            A = A_temp;
            B = B_temp;
            C = C_temp;
            D = D_temp;
            G = G_temp;
            alpha = alpha_temp;
            x0 = x0_temp;
        end
    end
    
    % Return cell arrays if multiple networks requested
    if num_networks > 1
        A = A_all;
        B = B_all;
        C = C_all;
        D = D_all;
        G = G_all;
        alpha = alpha_all;
        x0 = x0_all;
    end
end

function is_similar = isSimilarNetwork(A1, B1, C1, G1, A2, B2, C2, G2)
    % Check if two networks are too similar
    % Returns true if networks are nearly identical
    
    threshold = 0.01;  % Similarity threshold
    
    % Check if matrices have same sparsity pattern and similar values
    A_diff = norm(A1 - A2, 'fro') / (norm(A1, 'fro') + eps);
    B_diff = norm(B1 - B2, 'fro') / (norm(B1, 'fro') + eps);
    C_diff = norm(C1 - C2, 'fro') / (norm(C1, 'fro') + eps);
    G_diff = norm(G1 - G2, 'fro') / (norm(G1, 'fro') + eps);
    
    % Networks are similar if all differences are below threshold
    is_similar = (A_diff < threshold) && (B_diff < threshold) && ...
                 (C_diff < threshold) && (G_diff < threshold);
end

function states_distribution = distribute_states(n, num_subnetworks, min_states, max_states)
    % Distribute n states among num_subnetworks
    % Each subnetwork gets between min_states and max_states
    
    states_distribution = zeros(num_subnetworks, 1);
    remaining_states = n;
    
    for i = 1:num_subnetworks-1
        % Calculate how many states this subnet can have
        min_allowed = min_states;
        max_allowed = min(max_states, remaining_states - (num_subnetworks - i) * min_states);
        
        if max_allowed < min_allowed
            max_allowed = min_allowed;
        end
        
        % Assign random number of states within bounds
        states_distribution(i) = randi([min_allowed, max_allowed]);
        remaining_states = remaining_states - states_distribution(i);
    end
    
    % Last subnetwork gets remaining states
    states_distribution(num_subnetworks) = remaining_states;
    
    % Ensure last subnet has at least 1 state
    if states_distribution(num_subnetworks) < 1
        % Steal from previous subnetworks
        for i = (num_subnetworks-1):-1:1
            if states_distribution(i) > min_states
                transfer = min(states_distribution(i) - min_states, 1 - states_distribution(num_subnetworks));
                states_distribution(i) = states_distribution(i) - transfer;
                states_distribution(num_subnetworks) = states_distribution(num_subnetworks) + transfer;
                if states_distribution(num_subnetworks) >= 1
                    break;
                end
            end
        end
    end
end

function reaches_input = checkReachability(A, B, G, start_state)
    % Quick backward reachability check
    reachable_states = start_state;
    to_explore = start_state;
    explored = [];
    
    max_iterations = size(A, 1);
    iteration = 0;
    
    while ~isempty(to_explore) && iteration < max_iterations
        iteration = iteration + 1;
        current = to_explore(1);
        to_explore(1) = [];
        
        if ismember(current, explored)
            continue;
        end
        explored = [explored, current];
        
        % Find predecessors
        predecessors = find(A(:, current) ~= 0);
        
        for pred = predecessors'
            if ~ismember(pred, reachable_states)
                reachable_states = [reachable_states, pred];
                to_explore = [to_explore, pred];
            end
        end
        
        % Check bidirectional
        successors = find(A(current, :) ~= 0);
        for succ = successors
            if A(succ, current) ~= 0 && ~ismember(succ, reachable_states)
                reachable_states = [reachable_states, succ];
                to_explore = [to_explore, succ];
            end
        end
    end
    
    % Check if any input reaches these states
    reaches_input = false;
    for state = reachable_states
        if any(B(state, :) ~= 0) || any(G(state, :) ~= 0)
            reaches_input = true;
            break;
        end
    end
end