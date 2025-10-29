%% FON_Dtilde
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