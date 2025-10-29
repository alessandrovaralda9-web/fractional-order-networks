%% FON_z_sim
% Compute the convolution of psi functions with input states
% Inputs:
%   x     - Matrix of states [n,N] (n states, N time steps)
%   alpha - Vector of alpha values [alpha_1, alpha_2, ...]
%   J     - Number of terms to include in the convolution sum
% Outputs:
%   z     - Matrix of computed z values [n,N-J+1]
%           Note: z is computed only for time steps where all required x values are available

function z = FON_z_sim(x, alpha, J)

% Validate inputs
validateattributes(x, {'numeric'}, {'2d'});
validateattributes(alpha, {'numeric'}, {'vector'});
validateattributes(J, {'numeric'}, {'scalar', 'positive', 'integer'});

[n, N] = size(x);

% Check if alpha length matches the number of states
if length(alpha) ~= n
    error('The length of alpha must match the number of states n');
end

% Check if we have enough data points
if N < J
    error('Not enough time steps in x to compute z with convolution length J');
end

% Initialize output (we can only compute N-J+1 values of z)
z = zeros(n, N-J+1);

% Compute D matrices for j=0 to J-1
% Note: D{1} corresponds to j=0, D{2} to j=1, etc.
% For j=0, D is the identity matrix (not computed by FON_Dtilde)
D = cell(J, 1);
D{1} = eye(n); % D^(0,alpha)_0 is the identity matrix

% Compute the remaining D matrices using FON_Dtilde
if J > 1
    % Use FON_Dtilde with zero delta_alpha to compute D^(0,alpha)_j for j>0
    % We're using a trick: D^(0,alpha)_j is equivalent to Dtilde(0,alpha,j) - Dtilde(0,0,j)
    % Since Dtilde(0,0,j) is zero, we just need Dtilde(0,alpha,j)
    D_from_func = FON_Dtilde(zeros(size(alpha)), alpha, J-1);

    % Copy the computed D matrices (starting from j=1)
    for j = 1:J-1
        D{j+1} = D_from_func{j};
    end
end

for k = 0 : N-2 % Z[k]

    max_id = min(k+1,J-1);
    z(:,k+1) = zeros(size(x(:,k+1)));
    D_matrix = [];
    X_vect = [];
    for i = 1 : max_id+1
        D_matrix = [D_matrix , D{i}];
        X_vect = [X_vect ; x(:,max_id+2-i)];
    end
    z(:,k+1) = D_matrix*X_vect;

end

end
