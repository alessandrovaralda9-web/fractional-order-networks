function [bestAlpha, bestA, estError] = no_input_estimation_with_known_params(x, tau, A_known, A2est, alpha_known, alpha2est)
    % NO_INPUT_ESTIMATION_WITH_KNOWN_PARAMS - Estimate A and alpha with partial knowledge
    %
    % This function performs gradient descent optimization with partially known 
    % A matrix and alpha vector, reusing known values to reduce computation.
    %
    % Inputs:
    %   x           - State trajectory data (n × T)
    %   tau         - Memory parameter for FON_z_sim
    %   A_known     - Matrix with known values of A (n × n)
    %   A2est       - Binary mask: 1 = estimate, 0 = use known (n × n)
    %   alpha_known - Vector with known alpha values (n × 1) [OPTIONAL]
    %   alpha2est   - Binary mask: 1 = estimate alpha, 0 = use known (n × 1) [OPTIONAL]
    %
    % Outputs:
    %   bestAlpha   - Estimated fractional orders (n × 1)
    %   bestA       - Estimated state matrix (n × n)
    %   estError    - Frobenius norm of estimation error
    
    n = size(x, 1);
    
    % Handle optional alpha parameters
    if nargin < 5 || isempty(alpha_known)
        alpha_known = 0.5 * ones(n, 1);
        alpha2est = ones(n, 1);  % Estimate all by default
    end
    
    if nargin < 6 || isempty(alpha2est)
        alpha2est = ones(n, 1);  % Estimate all by default
    end
    
    % Find which alpha values need to be estimated
    unknownAlphaIdx = find(alpha2est == 1);
    
    % Special case: all alphas are known, just solve for A
    if isempty(unknownAlphaIdx)
        bestAlpha = alpha_known;
        z = FON_z_sim(x, bestAlpha, tau);
        k = size(x, 2) - 1;
        X = x(:, 1:k);
        Z = z(:, 1:k);
        bestA = solveForA_withKnown(Z, X, A_known, A2est);
        estZ = bestA * X;
        estError = norm(Z - estZ, 'fro');
        return;
    end
    
    % Optimization setup - only optimize unknown alphas
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off');
    
    % Initial guess for unknown alphas (use provided known values as starting point)
    initialAlphaUnknown = alpha_known(unknownAlphaIdx);
    
    % Objective function that only varies unknown alphas
    objectiveFunction = @(alpha_unknown) computeModelError_masked(...
        alpha_unknown, x, tau, A_known, A2est, alpha_known, unknownAlphaIdx);
    
    % Perform optimization over unknown alphas only
    [bestAlphaUnknown, ~, exit_flag] = fminunc(objectiveFunction, initialAlphaUnknown, options);
    
    % Reconstruct full alpha vector
    bestAlpha = alpha_known;
    bestAlpha(unknownAlphaIdx) = bestAlphaUnknown;
    
    % Compute final A and errors using the best alpha found
    z = FON_z_sim(x, bestAlpha, tau);
    k = size(x, 2) - 1;
    X = x(:, 1:k);
    Z = z(:, 1:k);
    bestA = solveForA_withKnown(Z, X, A_known, A2est);
    estZ = bestA * X;
    estError = norm(Z - estZ, 'fro');
end

function error = computeModelError_masked(alpha_unknown, x, tau, A_known, A2est, alpha_known, unknownAlphaIdx)
    % COMPUTEMODELERROR_MASKED - Compute error for given unknown alphas
    %
    % This function evaluates the reconstruction error when only some alpha
    % values are being optimized, with others held fixed.
    
    % Reconstruct full alpha vector from known and unknown parts
    alpha = alpha_known;
    alpha(unknownAlphaIdx) = alpha_unknown;
    
    % Compute Z using full alpha vector
    z = FON_z_sim(x, alpha, tau);
    k = size(x, 2) - 1;
    X = x(:, 1:k);
    Z = z(:, 1:k);
    
    % Solve for A incorporating known values
    A = solveForA_withKnown(Z, X, A_known, A2est);
    
    % Calculate error
    estZ = A * X;
    error = norm(Z - estZ, 'fro');
end



function A = solveForA_withKnown(Z, X, A_known, A2est)
    % SOLVEFORA_WITHKNOWN - Solve for A matrix with partially known entries
    %
    % This function efficiently solves for unknown entries of A by:
    %   1. Moving known contributions to the left side (adjusting Z)
    %   2. Solving only for unknown parameters (reducing problem size)
    %
    % Inputs:
    %   Z        - Target matrix (fractional derivative of states) (n × T)
    %   X        - State matrix (n × T)
    %   A_known  - Known values of A (n × n)
    %   A2est    - Binary mask: 1 = estimate, 0 = use known (n × n)
    %
    % Output:
    %   A        - Estimated A matrix (n × n)
    
    [m, p] = size(A_known);
    A = A_known;  % Start with known values
    
    for i = 1:m
        % Find unknown elements in this row (where A2est == 1)
        unknownColumns = find(A2est(i, :) == 1);
        
        if isempty(unknownColumns)
            % All values known for this row, nothing to estimate
            continue;
        end
        
        % Find known elements (where A2est == 0)
        knownColumns = find(A2est(i, :) == 0);
        
        % Adjust Z by subtracting the contribution from known elements
        % This moves the known part to the left side: Z_adj = Z - A_known*X
        if ~isempty(knownColumns)
            Z_adjusted = Z(i, :) - A_known(i, knownColumns) * X(knownColumns, :);
        else
            % All elements are unknown
            Z_adjusted = Z(i, :);
        end
        
        % Solve for unknown elements only (reduced problem size)
        X_unknown = X(unknownColumns, :);
        
        % Check rank and solve
        if rank(X_unknown) == size(X_unknown, 1)
            a_row_unknown = Z_adjusted / X_unknown;
        else
            a_row_unknown = Z_adjusted * pinv(X_unknown);
        end
        
        % Place estimated values in A
        A(i, unknownColumns) = a_row_unknown;
    end
end