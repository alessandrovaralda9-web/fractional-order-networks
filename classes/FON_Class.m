classdef FON_Class < handle
    % FON_Class Fractional Order Network implementation
    %   This class implements a Fractional Order Network system with methods for
    %   simulation and analysis. It handles state-space representations with
    %   fractional order dynamics.

    properties (Access = public)
        % System matrices
        A           % State matrix
        B           % Input matrix
        C           % Output matrix
        D           % Feedthrough matrix

        % System parameters
        alpha = []  % Fractional order
        n = 0       % Number of states
        m = 1       % Number of inputs
        Ts = 1      % Sampling time

        % Simulation results
        t = []      % Time vector
        x = []      % State trajectories
        z = []      % Fractional state trajectories
        y = []      % Output trajectories

        % System properties
        stable = false      % Stability flag
        eig_T = []         % System eigenvalues
    end

    methods

        function obj = FON_Class(ss1, alpha, Ts)
            % FON_Class Constructor
            % Inputs:
            %   ss1   - State space system
            %   alpha - Fractional order (optional)
            %   Ts    - Sampling time (optional)
            % Outputs:
            %   obj   - FON_Class object initialized with the specified parameters
            %          The object contains:
            %          - System matrices (A, B, C, D)
            %          - Fractional order (alpha)
            %          - System dimensions (n, m)
            %          - Sampling time (Ts)
            %          - Stability information (stable, eig_T)

            validateattributes(ss1, {'ss'}, {'nonempty'}, 'FON_Class', 'ss1', 1)

            % Initialize system matrices
            obj.A = ss1.A;
            obj.B = ss1.B;
            obj.C = ss1.C;
            obj.D = ss1.D;
            obj.n = size(ss1.A, 1);

            if nargin >= 2
                validateattributes(alpha, {'numeric'}, {'vector'}, 'FON_Class', 'alpha', 2)
                obj.alpha = alpha;
                % Get first D matrix from cell array for stability check
                D_cells = obj.FON_Dtilde(zeros(1, obj.n), alpha, 1);
                obj.eig_T = eig(ss1.A - D_cells{1});
                obj.stable = all(abs(obj.eig_T) < 1);
            end

            if nargin >= 3
                validateattributes(Ts, {'numeric'}, {'positive', 'scalar'}, 'FON_Class', 'Ts', 3)
                obj.Ts = Ts;
            end
        end

        % checked
        % Simulate the FON system dynamics
        % Inputs:
        %   u  - Input sequence
        %   x0 - Initial state
        %   v  - State noise parameters (struct with mu and sigma)
        %   w  - Output noise parameters (struct with mu and sigma)
        %   J  - Memory length parameter
        % Outputs:
        %   Updates obj.x - State trajectories matrix (n × N)
        %   Updates obj.y - Output trajectories matrix (n × (N-1))
        function fsim(obj, u, x0, v, w, J)

            validateattributes(u, {'numeric'}, {'2d'}, 'fsim', 'u')
            validateattributes(x0, {'numeric'}, {'vector', 'numel', obj.n}, 'fsim', 'x0')
            validateattributes(J, {'numeric'}, {'positive', 'integer', 'scalar'}, 'fsim', 'J')

            % Initialize simulation arrays
            N = length(u(1,:));
            obj.x = zeros(obj.n, N);
            obj.y = zeros(size(obj.C, 1), N-1);

            % Set initial state
            obj.x(:,1) = x0;

            % Pre-compute D matrices for the entire simulation
            D_cells = obj.FON_Dtilde(zeros(size(obj.alpha)), obj.alpha, min(N,J));

            % Simulation loop
            for i = 2:N+1
                max_el = min(i-1, round(J));

                % Calculate fractional sum term
                frac_sum = zeros(obj.n, 1);
                for j = 1:max_el
                    if i-j > 0
                        frac_sum = frac_sum + D_cells{j} * obj.x(:,i-j);
                    end
                end

                % State update with noise
                obj.x(:,i) = obj.A * obj.x(:,i-1) + obj.B * u(:,i-1) - frac_sum + ...
                    normrnd(v.mu, v.sigma, size(obj.x(:,1)));

                % Output equation with noise
                obj.y(:,i-1) = obj.C * obj.x(:,i-1) + obj.D * u(:,i-1) + ...
                    normrnd(w.mu, w.sigma, size(obj.C * obj.x(:,i-1)));
            end
        end

        % checked
        % Compute input sequence for shifting alpha value
        % Inputs:
        %   u           - Original input sequence
        %   x0          - Initial state
        %   v           - State noise parameters
        %   w           - Output noise parameters
        %   alpha_des   - Desired alpha value
        %   shift_index - Index to start shift
        % Outputs:
        %   u_tilde    - Modified input sequence (n × N)
        %   Updates obj.x - State trajectories matrix (n × N)
        %   Updates obj.y - Output trajectories matrix (n × N)
        function u_tilde = fsim_shift(obj, u, x0, v, w, alpha_des, shift_index, J)

            %validateattributes(alpha_des, {'numeric'}, {'vector'}, 'fsim_shift', 'alpha_des')
            %validateattributes(shift_index, {'numeric'}, {'positive', 'integer'}, 'fsim_shift', 'shift_index')

            validateattributes(u, {'numeric'}, {'2d'}, 'fsim', 'u')
            validateattributes(x0, {'numeric'}, {'vector', 'numel', obj.n}, 'fsim', 'x0')
            validateattributes(J, {'numeric'}, {'positive', 'integer', 'scalar'}, 'fsim', 'J')

            % Initialize arrays
            N = length(u(1,:));
            obj.x = zeros(obj.n, N);
            obj.y = zeros(size(obj.C, 1), N);
            u_tilde = zeros(size(u));

            % Set initial state
            obj.x(:,1) = x0;

            % Compute alpha difference
            delta_alpha = alpha_des - obj.alpha;

            % Pre-compute D matrices for both current and shift cases
            D_cells = obj.FON_Dtilde(zeros(size(obj.alpha)), obj.alpha, N);
            D_cells_shift = obj.FON_Dtilde(obj.alpha, delta_alpha, N);

            % Simulation loop
            for i = 2:N+1

                max_el = min(i-1, round(J));
                % Calculate current and shift fractional sums
                frac_sum = zeros(obj.n, 1);
                shift_sum = zeros(obj.n, 1);

                for j = 1:max_el
                    if i-j > 0
                        frac_sum = frac_sum + D_cells{j} * obj.x(:,i-j);
                        shift_sum = shift_sum + D_cells_shift{j} * obj.x(:,i-j);
                    end
                end

                % Compute input
                if i-1 < shift_index
                    u_tilde(:,i-1) = u(:,i-1);
                else
                    u_tilde(:,i-1) = u(:,i-1) - obj.B \ shift_sum;
                end

                % State update with noise
                obj.x(:,i) = obj.A * obj.x(:,i-1) + obj.B * u_tilde(:,i-1) - frac_sum + ...
                    normrnd(v.mu, v.sigma, size(obj.x(:,1)));

                % Output equation with noise
                obj.y(:,i) = obj.C * obj.x(:,i-1) + obj.D * u_tilde(:,i-1) + ...
                    normrnd(w.mu, w.sigma, size(obj.C * obj.x(:,i-1)));
            end
        end

        % checked
        % Simulate the FON system dynamics
        % Inputs:
        %   u  - Input sequence
        %   x0 - Initial state
        %   v  - State noise parameters (struct with mu and sigma)
        %   w  - Output noise parameters (struct with mu and sigma)
        %   J  - Memory length parameter
        % Outputs:
        %   Updates obj.x - State trajectories matrix (n × N)
        %   Updates obj.y - Output trajectories matrix (n × (N-1))
        function fsim_statefeedback(obj, K, x0, v, w, J, N)

            validateattributes(x0, {'numeric'}, {'vector', 'numel', obj.n}, 'fsim', 'x0')
            validateattributes(J, {'numeric'}, {'positive', 'integer', 'scalar'}, 'fsim', 'J')

            % Initialize simulation arrays
            obj.x = zeros(obj.n, N);
            obj.y = zeros(size(obj.C, 1), N-1);

            % Set initial state
            obj.x(:,1) = x0;

            % Pre-compute D matrices for the entire simulation
            D_cells = obj.FON_Dtilde(zeros(size(obj.alpha)), obj.alpha, min(N,J));

            % Simulation loop
            for i = 2:N+1
                max_el = min(i-1, round(J));

                % Calculate fractional sum term
                frac_sum = zeros(obj.n, 1);
                for j = 1:max_el
                    if i-j > 0
                        frac_sum = frac_sum + D_cells{j} * obj.x(:,i-j);
                    end
                end

                % State update with noise
                obj.x(:,i) = (obj.A-obj.B*K) * obj.x(:,i-1) - frac_sum + ...
                    normrnd(v.mu, v.sigma, size(obj.x(:,1)));

                % Output equation with noise
                obj.y(:,i-1) = (obj.C-obj.D*K) * obj.x(:,i-1) + ...
                    normrnd(w.mu, w.sigma, size(obj.C * obj.x(:,i-1)));
            end
        end

        % checked
        % fsim_servo Simulate the FON system with servo control
        %
        % This function implements servo control for fractional-order systems where:
        % u[k] = Lr*r[k] - K*x[k]
        %
        % Inputs:
        %   K   - State feedback gain matrix
        %   Lr  - Reference feedforward filter matrix
        %   r   - Reference signal sequence (m × N matrix or vector)
        %   x0  - Initial state vector
        %   v   - State noise parameters (struct with mu and sigma)
        %   w   - Output noise parameters (struct with mu and sigma)
        %   J   - Memory length parameter
        %
        % Outputs:
        %   Updates obj.x - State trajectories matrix (n × N)
        %   Updates obj.y - Output trajectories matrix (p × (N-1))
        %   Updates obj.u - Control input trajectories matrix (m × (N-1))
        function u = fsim_servo(obj, K, Lr, r, x0, v, w, J)
            % fsim_servo Simulate the FON system with servo control
            %
            % This function implements servo control for fractional-order systems where:
            % u[k] = Lr*r[k] - K*x[k]
            %
            % Inputs:
            %   K   - State feedback gain matrix
            %   Lr  - Reference feedforward filter matrix
            %   r   - Reference signal sequence (m × N matrix or vector)
            %   x0  - Initial state vector
            %   v   - State noise parameters (struct with mu and sigma)
            %   w   - Output noise parameters (struct with mu and sigma)
            %   J   - Memory length parameter
            %
            % Outputs:
            %   u       - Control input trajectories matrix (m × (N-1))
            %   Updates obj.x - State trajectories matrix (n × N)
            %   Updates obj.y - Output trajectories matrix (p × (N-1))

            % Input validation
            validateattributes(x0, {'numeric'}, {'vector', 'numel', obj.n}, 'fsim_servo', 'x0');
            validateattributes(J, {'numeric'}, {'positive', 'integer', 'scalar'}, 'fsim_servo', 'J');
            validateattributes(K, {'numeric'}, {'2d'}, 'fsim_servo', 'K');
            validateattributes(Lr, {'numeric'}, {'2d'}, 'fsim_servo', 'Lr');
            validateattributes(r, {'numeric'}, {}, 'fsim_servo', 'r');

            % Handle reference signal dimensions
            if isvector(r)
                r = r(:)';  % Make it a row vector
                N = length(r);
                r = repmat(r, size(obj.B, 2), 1);  % Replicate for each input if needed
            else
                N = size(r, 2);
            end

            % Verify matrix dimensions
            if size(K, 2) ~= obj.n
                error('K matrix must have %d columns to match system states', obj.n);
            end
            if size(K, 1) ~= size(obj.B, 2)
                error('K matrix must have %d rows to match system inputs', size(obj.B, 2));
            end
            if size(Lr, 1) ~= size(obj.B, 2)
                error('Lr matrix must have %d rows to match system inputs', size(obj.B, 2));
            end
            if size(Lr, 2) ~= size(r, 1)
                error('Lr matrix columns must match reference signal dimensions');
            end

            % Initialize simulation arrays
            obj.x = zeros(obj.n, N);
            obj.y = zeros(size(obj.C, 1), N-1);

            % Initialize control input array
            u = zeros(size(obj.B, 2), N-1);

            % Set initial state
            obj.x(:,1) = x0;

            % Pre-compute D matrices for the entire simulation
            D_cells = obj.FON_Dtilde(zeros(size(obj.alpha)), obj.alpha, min(N,J));

            % Simulation loop
            for i = 2:N+1
                max_el = min(i-1, round(J));

                % Calculate fractional sum term
                frac_sum = zeros(obj.n, 1);
                for j = 1:max_el
                    if i-j > 0
                        frac_sum = frac_sum + D_cells{j} * obj.x(:,i-j);
                    end
                end

                % Compute control input: u[k] = Lr*r[k] - K*x[k]
                u_k = Lr * r(:,i-1) - K * obj.x(:,i-1);
                u(:,i-1) = u_k;

                % State update with servo control and noise
                % x[k+1] = A*x[k] + B*u[k] - frac_sum + noise
                % x[k+1] = A*x[k] + B*(Lr*r[k] - K*x[k]) - frac_sum + noise
                % x[k+1] = (A - B*K)*x[k] + B*Lr*r[k] - frac_sum + noise
                obj.x(:,i) = obj.A * obj.x(:,i-1) + obj.B * u_k - frac_sum + ...
                    normrnd(v.mu, v.sigma, size(obj.x(:,1)));

                % Output equation with servo control and noise
                % y[k] = C*x[k] + D*u[k] + noise
                % y[k] = C*x[k] + D*(Lr*r[k] - K*x[k]) + noise
                % y[k] = (C - D*K)*x[k] + D*Lr*r[k] + noise
                obj.y(:,i-1) = obj.C * obj.x(:,i-1) + obj.D * u_k + ...
                    normrnd(w.mu, w.sigma, size(obj.C * obj.x(:,i-1)));
            end
        end


        % checked
        % Compute input sequence for shifting alpha value WITH state feedback
        % Inputs:
        %   u           - Original input sequence
        %   K           - State feedback gain matrix
        %   x0          - Initial state
        %   v           - State noise parameters
        %   w           - Output noise parameters
        %   alpha_des   - Desired alpha value
        %   shift_index - Index to start shift
        %   J           - Memory length parameter
        % Outputs:
        %   u_tilde    - Modified input sequence (n × N)
        %   Updates obj.x - State trajectories matrix (n × N)
        %   Updates obj.y - Output trajectories matrix (n × N)
        function u_tilde = fsim_shift_statefeedback(obj, u, K, x0, v, w, alpha_des, shift_index, J)

            validateattributes(u, {'numeric'}, {'2d'}, 'fsim_shift_statefeedback', 'u')
            validateattributes(K, {'numeric'}, {'2d'}, 'fsim_shift_statefeedback', 'K')
            validateattributes(x0, {'numeric'}, {'vector', 'numel', obj.n}, 'fsim_shift_statefeedback', 'x0')
            validateattributes(J, {'numeric'}, {'positive', 'integer', 'scalar'}, 'fsim_shift_statefeedback', 'J')

            % Initialize arrays
            N = length(u(1,:));
            obj.x = zeros(obj.n, N);
            obj.y = zeros(size(obj.C, 1), N);
            u_tilde = zeros(size(u));

            % Set initial state
            obj.x(:,1) = x0;

            % Compute alpha difference
            delta_alpha = alpha_des - obj.alpha;

            % Pre-compute D matrices for both current and shift cases
            D_cells = obj.FON_Dtilde(zeros(size(obj.alpha)), obj.alpha, N);
            D_cells_shift = obj.FON_Dtilde(obj.alpha, delta_alpha, N);

            % Simulation loop
            for i = 2:N+1

                max_el = min(i-1, round(J));
                % Calculate current and shift fractional sums
                frac_sum = zeros(obj.n, 1);
                shift_sum = zeros(obj.n, 1);

                for j = 1:max_el
                    if i-j > 0
                        frac_sum = frac_sum + D_cells{j} * obj.x(:,i-j);
                        shift_sum = shift_sum + D_cells_shift{j} * obj.x(:,i-j);
                    end
                end

                % Compute combined input (feedforward + alpha shift + state feedback)
                if i-1 < shift_index
                    % Before shift index: only apply state feedback
                    u_tilde(:,i-1) = u(:,i-1) - K * obj.x(:,i-1);
                else
                    % After shift index: apply both alpha shift and state feedback
                    u_tilde(:,i-1) = u(:,i-1) - obj.B \ shift_sum - K * obj.x(:,i-1);
                end

                % State update with noise
                obj.x(:,i) = obj.A * obj.x(:,i-1) + obj.B * u_tilde(:,i-1) - frac_sum + ...
                    normrnd(v.mu, v.sigma, size(obj.x(:,1)));

                % Output equation with noise
                obj.y(:,i) = obj.C * obj.x(:,i-1) + obj.D * u_tilde(:,i-1) + ...
                    normrnd(w.mu, w.sigma, size(obj.C * obj.x(:,i-1)));
            end
        end



        function z_sim(obj, alpha, J_perc)
            % Simulate fractional behavior of the network
            % Inputs:
            %   alpha   - Fractional order
            %   J_perc  - Memory length percentage
            % Outputs:
            %   Updates obj.z - Fractional state trajectories matrix (n × (N-1))

            validateattributes(J_perc, {'numeric'}, {'positive', 'scalar', '<=', 100}, 'z_sim', 'J_perc')

            N = size(obj.x, 2) - 1;
            obj.z = zeros(size(obj.x, 1), N);
            J = J_perc/100;

            % Pre-compute psi values for each alpha_i
            psi_vect = zeros(N+1, obj.n);
            for i = 1:obj.n
                psi_vect(:,i) = FON_Class.FON_psi(alpha(i), N);
            end

            for i = 1:N
                l = min([i, round(N*J)]);
                z_sum = zeros(size(obj.x, 1), 1);

                for j = 1:l
                    z_sum = z_sum + psi_vect(j+1,:)' .* obj.x(:,l+1-j);
                end

                obj.z(:,i) = obj.x(:,i+1) + z_sum;
            end
        end


        function G = FON_Gmatrix(obj, alpha_sequence, k1, q1)
            % Computes the state transition matrix G_k1^(0,α[q1]) for Fractional Order Networks
            % as defined in "Steering Fractional Exponents in Fractional-Order Networks..."
            %
            % Inputs:
            %   alpha_sequence  - Matrix where each row contains alpha vector for a time step
            %                     alpha_sequence(i,:) contains α[i-1]
            %   k1              - Time index k1 (equivalent to k+1 in the paper)
            %   q1              - Index q1 (equivalent to q+1 in the paper)
            %
            % Output:
            %   G               - State transition matrix G_k1^(0,α[q1])

            % Base case: k1 = 0
            if k1 == 0
                G = eye(obj.n);
                return;
            end

            % Get the relevant alpha vector
            alpha_q1 = alpha_sequence(q1+1,:)';

            % Initialize result
            G = zeros(obj.n);

            % Pre-compute D_tilde matrices for efficiency
            D_cells = obj.FON_Dtilde(zeros(size(alpha_q1)), alpha_q1, k1);

            % Compute recursively according to the definition
            for j = 0:k1-1
                % Determine A_tilde based on j value (matches paper definition)
                if j == 0
                    A_tilde = obj.A - D_cells{j+1};
                else
                    A_tilde = -D_cells{j+1};
                end

                % Recursively compute G with proper indices
                if q1-j >= 0
                    G_prev = obj.FON_Gmatrix(alpha_sequence, k1-j-1, q1-j);
                    % Add to sum
                    G = G + A_tilde * G_prev;
                end
            end
        end


        function C = FON_Cntrl(obj, alpha_sequence, k1)
            % FON_Cntrl Computes the fractional controllability matrix C^(0,α[k1])_k1-1
            % as defined in "Steering Fractional Exponents in Fractional-Order Networks..."
            %
            % Inputs:
            %   alpha_sequence - Matrix where each row contains alpha vector for a time step
            %                    alpha_sequence(i,:) contains α[i-1]
            %   k1             - Time index k1 (equivalent to k+1 in the paper)
            %
            % Output:
            %   C              - Fractional controllability matrix C^(0,α[k1])_(k1-1)

            % Verify inputs
            if k1 < 1
                error('k1 must be at least 1');
            end

            if size(alpha_sequence, 1) < k1
                error('alpha_sequence must have at least k1 rows');
            end

            % Get dimensions
            n = obj.n;
            m = size(obj.B, 2);

            % Initialize the controllability matrix
            C = zeros(n, k1*m);

            % Fill the controllability matrix according to the definition
            % C^(0,α[k1])_(k1-1) = [G^(0,α[k1])_0 B  G^(0,α[k1])_1 B  ...  G^(0,α[k1])_(k1-1) B]
            for i = 0:k1-1
                % Compute the G matrix for the current index
                G = obj.FON_Gmatrix(alpha_sequence, i, k1);

                % Place G*B in the appropriate position in the controllability matrix
                C(:, i*m+1:(i+1)*m) = G * obj.B;
            end
        end


        function G_cells = FON_Gmatrix_all(obj, alpha_d, k)
            % FON_Gmatrix_all Computes all state transition matrices G_0, G_1, ..., G_k, G_{k+1}
            % for constant fractional order alpha_d using efficient matrix formulation
            %
            % Uses: G_{k+1} = [A_tilde_1, A_tilde_2, ..., A_tilde_{k+1}] * [G_k; G_{k-1}; ...; G_0]
            %
            % Inputs:
            % alpha_d - Constant fractional order vector (n_alpha x 1)
            % k       - Maximum index (will compute up to G_{k+1})
            %
            % Output:
            % G_cells - Cell array where G_cells{i} = G_{i-1}^(0,α_d)
            %          G_cells{1} = G_0, G_cells{2} = G_1, ..., G_cells{k+2} = G_{k+1}

            % Validate inputs
            validateattributes(alpha_d, {'numeric'}, {'vector'}, 'FON_Gmatrix_all', 'alpha_d')
            validateattributes(k, {'numeric'}, {'nonnegative', 'integer', 'scalar'}, 'FON_Gmatrix_all', 'k')

            % Ensure column vector
            alpha_d = alpha_d(:);

            % Initialize cell array
            G_cells = cell(k+2, 1);

            % Base case: G_0 = I
            G_cells{1} = eye(obj.n);

            % If k < 0, only return G_0
            if k < 0
                G_cells = {eye(obj.n)};
                return;
            end

            % Pre-compute all needed D_tilde matrices
            D_cells = obj.FON_Dtilde(zeros(size(alpha_d)), alpha_d, k+1);

            % Initialize G_stack with G_0 (column-stacked format)
            G_stack = G_cells{1};  % (n x n)

            % Iteratively compute G_1, G_2, ..., G_{k+1}
            for i = 1:k+1  % Computing G_i where i goes from 1 to k+1

                % Build A_tilde matrix: [A_tilde_1, A_tilde_2, ..., A_tilde_i]
                A_tilde_row = [];

                for j = 1:i
                    if j == 1
                        A_tilde_j = obj.A - D_cells{j};  % A - D_1
                    else
                        A_tilde_j = -D_cells{j};         % -D_j
                    end
                    A_tilde_row = [A_tilde_row, A_tilde_j];  % Horizontal concatenation
                end

                % Compute G_i = A_tilde_row * G_stack
                % where G_stack = [G_{i-1}; G_{i-2}; ...; G_0] (column-stacked)
                G_cells{i+1} = A_tilde_row * G_stack;

                % Update G_stack: prepend G_i to get [G_i; G_{i-1}; ...; G_0]
                G_stack = [G_cells{i+1}; G_stack];  % Vertical concatenation
            end
        end

        function u_sequence = FON_OptimalInput(obj, x0, x_desired, alpha_desired, k)
            % FON_OptimalInput Computes the optimal input sequence to reach desired state
            % using the fractional controllability approach
            %
            % Inputs:
            %   x0             - Initial state vector (n x 1)
            %   x_desired      - Desired final state x[k+1] (n x 1)
            %   alpha_desired  - Desired fractional order vector (n_alpha x 1)
            %   k              - Time horizon (final time index)
            %
            % Output:
            %   u_sequence     - Optimal input sequence matrix (m x (k+1))
            %                    u_sequence(:,1) = u[k], u_sequence(:,2) = u[k-1], ..., u_sequence(:,k+1) = u[0]

            % Validate inputs
            validateattributes(x0, {'numeric'}, {'vector', 'numel', obj.n}, 'FON_OptimalInput', 'x0')
            validateattributes(x_desired, {'numeric'}, {'vector', 'numel', obj.n}, 'FON_OptimalInput', 'x_desired')
            validateattributes(alpha_desired, {'numeric'}, {'vector'}, 'FON_OptimalInput', 'alpha_desired')
            validateattributes(k, {'numeric'}, {'positive', 'integer', 'scalar'}, 'FON_OptimalInput', 'k')

            % Ensure column vectors
            x0 = x0(:);
            x_desired = x_desired(:);
            alpha_desired = alpha_desired(:);

            % Get dimensions
            n = obj.n;
            m = size(obj.B, 2);

            % OPTIMIZED APPROACH: Use FON_Gmatrix_all to get all G matrices at once
            % This is much more efficient than the recursive approach
            G_cells = obj.FON_Gmatrix_all(alpha_desired, k);

            % Construct controllability matrix C^(0,α_d)_k directly
            % C = [G_0*B | G_1*B | ... | G_k*B]
            C = zeros(n, (k+1)*m);
            for i = 0:k
                C(:, i*m+1:(i+1)*m) = G_cells{i+1} * obj.B;  % G_i is stored in G_cells{i+1}
            end

            % Extract G_{k+1}^(0,α_d) directly from the cell array
            G = G_cells{k+2};  % G_{k+1} is stored in G_cells{k+2}

            % Compute the target vector: x[k+1] - G_{k+1}^(0,α[k]) * x0
            target_vector = x_desired - G * x0;

            % Compute the pseudo-inverse of the controllability matrix
            C_pinv = pinv(C);

            % Compute the optimal input sequence vector
            % [u[k]; u[k-1]; ...; u[0]] = pinv(C) * target_vector
            u_vector = C_pinv * target_vector;

            % Reshape the result into matrix form (m x (k+1))
            % u_sequence(:,1) = u[k], u_sequence(:,2) = u[k-1], ..., u_sequence(:,k+1) = u[0]
            u_sequence = fliplr(reshape(u_vector, m, k+1));

        end


        % function [u_tilde, u] = FON_ComputeInput(obj, x0, xd, alpha_sequence)
        %     % FON_ComputeInput Computes input sequences to shift from (x0,alpha0) to (xd,alphad)
        %     % Based on Theorem 2 from "Steering Fractional Exponents in Fractional-Order Networks..."
        %     %
        %     % Inputs:
        %     %   x0             - Initial state vector
        %     %   xd             - Desired final state vector
        %     %   alpha_sequence - Matrix where each row contains alpha vector for a time step
        %     %                    alpha_sequence(i,:) contains α[i-1]
        %     %
        %     % Outputs:
        %     %   u_tilde        - Input sequence for defractionalized system (Eq. 13)
        %     %   u              - Actual input sequence for original system (Eq. 10)
        % 
        %     % Determine transition time T from alpha_sequence length
        %     T = size(alpha_sequence, 1) - 1;
        % 
        %     % Verify input dimensions
        %     n = obj.n;
        %     m = size(obj.B, 2);
        % 
        %     validateattributes(x0, {'numeric'}, {'vector', 'numel', n}, 'FON_ComputeInput', 'x0');
        %     validateattributes(xd, {'numeric'}, {'vector', 'numel', n}, 'FON_ComputeInput', 'xd');
        % 
        %     % Ensure column vectors
        %     x0 = x0(:);
        %     xd = xd(:);
        % 
        %     % Step 1: Compute fractional controllability matrix C^(0,αd)_(T-1)
        %     C = obj.FON_Cntrl(alpha_sequence, T);
        % 
        %     % Step 2: Compute G^(0,αd)_T
        %     % Note: Using T+1 for q1 to get α[T]
        %     G = obj.FON_Gmatrix(alpha_sequence, T, T);
        % 
        %     % Step 3: Check controllability
        %     rank_C = rank(C);
        %     if rank_C < n
        %         warning('System is not (α,x)-controllable. Rank of controllability matrix: %d/%d', rank_C, n);
        %     end
        % 
        %     % Step 4: Compute input sequence u_tilde using Eq. 13
        %     % u_tilde = C^(0,αd)_(T-1)† (xd - G^(0,αd)_T x0)
        % 
        %     % Compute pseudo-inverse of C
        %     C_pinv = pinv(C);
        % 
        %     % Compute desired state change
        %     delta_x = xd - G*x0;
        % 
        %     % Compute input sequence
        %     u_tilde_vec = C_pinv * delta_x;
        % 
        %     % Reshape into time sequence
        %     u_tilde = reshape(u_tilde_vec, m, T)';
        % 
        %     % Step 5: Compute actual input sequence u using Eq. 10
        %     u = zeros(size(u_tilde));
        % 
        %     % Define zero noise for simulation (deterministic case)
        %     v.mu = zeros(n, 1);
        %     v.sigma = zeros(n, 1);
        %     w.mu = zeros(size(obj.C, 1), 1);
        %     w.sigma = zeros(size(obj.C, 1), 1);
        %     J = T; % Use full memory length for accurate simulation
        % 
        %     % Simulate system evolution with u_tilde to get state trajectory
        %     obj.fsim(u_tilde', x0, v, w, J);
        %     x_sim = obj.x; % Store state trajectory for computing u
        % 
        %     for k = 1:T
        %         % Calculate fractional sum term for current time step
        %         frac_sum = zeros(n, 1);
        % 
        %         % Pre-compute all D matrices at once
        %         D_cells_k = obj.FON_Dtilde(zeros(size(alpha_sequence(1,:))), alpha_sequence(k,:), T);
        %         D_cells_k1 = obj.FON_Dtilde(zeros(size(alpha_sequence(1,:))), alpha_sequence(k+1,:), T);
        % 
        %         % The sum goes up to min(k,J) where J is memory length
        %         % We use min(k-1,T) since we're iterating from 1 to k-1
        %         for j = 1:k
        %             frac_sum = frac_sum + (D_cells_k1{j}-D_cells_k{j}) * x_sim(:,k-j+1);
        %         end
        % 
        %         % Compute actual input using Eq. 10: Bu[k] = Bu_tilde[k] - sum(D_tilde*x)
        %         % Note the correct sign is MINUS
        %         if m == n && rank(obj.B) == n
        %             % B is square and invertible
        %             u(k,:) = u_tilde(k,:) - (obj.B \ frac_sum)';
        %         else
        %             % B is not square or not invertible, use pseudo-inverse
        %             u(k,:) = u_tilde(k,:) - (pinv(obj.B) * frac_sum)';
        %         end
        %     end
        % end
        % 

        % checked
        % omega_matrix Generate omega matrix O of size (J×n) by (J×n)
        %
        % This function creates a block matrix where:
        % - Diagonal blocks: A - Delta_matrix(1)
        % - Super-diagonal blocks: -Delta_matrix(2), -Delta_matrix(3), etc.
        % - Sub-diagonal blocks: zeros
        %
        % Inputs:
        %   J - Memory length parameter (positive integer)
        %
        % Outputs:
        %   O - Omega matrix of size (J×n) by (J×n) where n is the number of states
        %
        % Matrix structure:
        % O = [A-Δ₁,  -Δ₂,  -Δ₃, ... ]
        %     [0,    A-Δ₁,  -Δ₂, ... ]
        %     [0,     0,   A-Δ₁, ... ]
        %     [⋮,     ⋮,     ⋮,   ⋱  ]
        function O = omega_matrix(obj, J)


            % Input validation
            validateattributes(J, {'numeric'}, {'positive', 'integer', 'scalar'}, 'omega_matrix', 'J');

            % Get system dimensions
            n = obj.n;

            % Pre-compute all Delta matrices using FON_Dtilde
            % This generates Delta_matrix(1), Delta_matrix(2), ..., Delta_matrix(J)
            D_cells = obj.FON_Dtilde(zeros(size(obj.alpha)), obj.alpha, J);

            % Initialize the omega matrix
            O = zeros(J*n, J*n);

            % Fill the omega matrix block by block
            for i = 1:J
                for j = i:J  % Only fill upper triangular part (including diagonal)
                    % Calculate the row and column indices for the current n×n block
                    row_idx = (i-1)*n + 1 : i*n;
                    col_idx = (j-1)*n + 1 : j*n;

                    if i == j
                        % Diagonal blocks: A - Delta_matrix(1)
                        O(row_idx, col_idx) = obj.A - D_cells{1};
                    else
                        % Super-diagonal blocks: -Delta_matrix(j-i+1)
                        delta_idx = j - i + 1;
                        if delta_idx <= J
                            O(row_idx, col_idx) = -D_cells{delta_idx};
                        end
                    end
                end
            end

            % Note: Sub-diagonal blocks remain zero (already initialized)
        end

    end



    methods (Static)

        % checked
        % Compute psi values for a scalar alpha
        % Inputs:
        %   alpha - Scalar fractional order
        %   j     - Number of iterations
        % Outputs:
        %   psi   - Column vector ((j+1) × 1) containing psi values
        %          starting from psi(alpha,0)
        function psi = FON_psi(alpha, j)

            validateattributes(alpha, {'numeric'}, {'scalar'});

            psi(1,1) = 1;
            for i = 1:j
                psi(i+1,1) = psi(i,1) * ((i-1-alpha)/i);
            end
        end

        % checked
        % Generate Dk matrices using the simplified formulation
        % Inputs:
        %   alpha0      - Vector of initial alpha values [alpha0_1, alpha0_2, ...]
        %   delta_alpha - Vector of alpha changes [delta_alpha_1, delta_alpha_2, ...]
        %   k           - Maximum time index
        % Outputs:
        %   D           - Cell array {k} where D{j} is the diagonal matrix D at time j
        function D = FON_Dtilde(alpha0, delta_alpha, k)

            validateattributes(alpha0, {'numeric'}, {'vector'});
            validateattributes(delta_alpha, {'numeric'}, {'vector', 'numel', length(alpha0)});

            n = length(alpha0);

            % Initialize storage for D matrices
            D = cell(k, 1);

            % Check if alpha0 is a zero vector
            if all(alpha0 == 0)
                % Only compute psi for alpha+delta_alpha
                psi_alpha_plus_delta = zeros(k+1, n);
                for i = 1:n
                    psi_alpha_plus_delta(:,i) = FON_Class.FON_psi(delta_alpha(i), k);
                end

                % For each time step j, compute D{j} = Dtilde(0,alpha+delta_alpha,j)
                for j = 1:k
                    D{j} = diag(psi_alpha_plus_delta(j+1,:));
                end
            else
                % Regular case: compute both terms
                psi_alpha = zeros(k+1, n);
                psi_alpha_plus_delta = zeros(k+1, n);

                % Compute psi for each alpha_i and (alpha_i + delta_alpha_i)
                for i = 1:n
                    psi_alpha(:,i) = FON_Class.FON_psi(alpha0(i), k);
                    psi_alpha_plus_delta(:,i) = FON_Class.FON_psi(alpha0(i) + delta_alpha(i), k);
                end

                % For each time step j, compute D{j} = Dtilde(0,alpha+delta_alpha,j) - Dtilde(0,alpha,j)
                for j = 1:k
                    D{j} = diag(psi_alpha_plus_delta(j+1,:) - psi_alpha(j+1,:));
                end
            end
        end
    end
end