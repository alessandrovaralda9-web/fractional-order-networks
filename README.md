# Fractional-Order Network Parameter Estimation with Hierarchical Decomposition

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)

## Overview

This repository contains the MATLAB implementation of parameter estimation methods for fractional-order networks (FONs) using hierarchical decomposition. The approach enables efficient identification of large-scale fractional-order systems by decomposing them into smaller, independently identifiable subsystems.

**Key Features:**
- Hierarchical decomposition of fractional-order networks based on graph structure
- Two-stage parameter estimation (fractional orders and system matrices)
- Validation framework with reconstruction from subsystem estimates
- Support for both controllable and uncontrollable inputs
- Comprehensive simulation and testing capabilities

This repository accompanies the paper:

> [Alessandro Varalda, Sergio Pequito]. "**Structural Identifiability in Fractional-Order Networks**" *[IFAC26]*, [2026]. [doi not yet available]

## Mathematical Framework

### System Model

The fractional-order network is described by:

```
Δ^α x[k+1] = A x[k] + B u[k] + G w[k]
y[k] = C x[k] + D u[k] + e[k]
```

Where:
- `x[k]` ∈ ℝⁿ: State vector
- `u[k]` ∈ ℝᵐ: Controllable inputs
- `w[k]` ∈ ℝᵛ: Uncontrollable inputs (disturbances)
- `y[k]` ∈ ℝᵒ: Output vector  
- `e[k]` ∈ ℝᵒ: Measurement noise
- `α` ∈ (0,1)ⁿ: Fractional orders for each state
- `Δ^α`: Grünwald-Letnikov fractional difference operator

### Parameter Estimation

**Stage 1: Zero-Input Estimation**
- Estimate `A` matrix and fractional orders `α` from zero-input response
- Uses gradient descent optimization with structural constraints

**Stage 2: Input-Driven Estimation**
- Estimate `B` and `G` matrices from input-output data
- Leverages estimated `A` and `α` from Stage 1

**Stage 3: Hierarchical Decomposition**
- Decompose network into subsystems based on input-output reachability
- Independently estimate parameters for each subsystem
- Reconstruct full system from subsystem estimates

## Requirements

### MATLAB Version
- MATLAB R2020a or later

### Required Toolboxes
- **Control System Toolbox** (for `ss` function)
- **Statistics and Machine Learning Toolbox** (for `corr` function)
- **Optimization Toolbox** (for `fminunc` function)

To check if you have these toolboxes:
```matlab
ver  % Lists all installed toolboxes
```

## Installation

1. Clone or download this repository:
   ```bash
   git clone https://github.com/yourusername/fractional-order-networks.git
   cd fractional-order-networks
   ```

2. Add the repository to your MATLAB path:
   ```matlab
   addpath(genpath('path/to/fractional-order-networks'))
   ```

---

## Repository structure

```
fractional-order-networks/
├── Fractional_Order_Network_param_est_with_subnet_decomp.m   % Main script
├── classes/
│   ├── FON_Class.m
│   └── FON_Graph_Class.m
├── functions/
│   ├── FON_Dtilde.m
│   ├── FON_psi.m
│   ├── FON_z_sim.m
│   ├── generateSparseNetwork.m
│   └── no_input_estimation_with_known_params.m
├── examples/
│   └── example_usage.m
├── README.md
├── LICENSE
└── .gitignore
```

> The file/folder list above matches the current repository contents. If you add or rename files later, please update this section.

---

## Quick Start

### Basic Usage

Run the main script to see the full parameter estimation workflow:

```matlab
% The main script performs 50 trials with random networks
run('Fractional_Order_Network_param_est_with_subnet_decomp.m')
```

This will:
1. Generate 50 random sparse fractional-order networks
2. Decompose each network hierarchically
3. Estimate parameters using full system (baseline)
4. Estimate parameters using subsystem decomposition
5. Compare and validate the approaches

### Simple Example

For a quick demonstration with a small network:

```matlab
% See examples/example_usage.m for a simplified demo
run('examples/example_usage.m')
```

### Custom Network

To test with your own network parameters:

```matlab
% Define system dimensions
n = 10;  % Number of states
m = 3;   % Number of controllable inputs
v = 2;   % Number of uncontrollable inputs
o = 3;   % Number of outputs

% Generate a single network
[A, B, C, D, G, alpha, x0] = generateSparseNetwork(n, m, v, o, 'high', 1);

% Create FON object
sys = ss(A, B, C, D);
fon = FON_Graph_Class(sys, alpha, G);

% Perform hierarchical decomposition
sub_systems = fon.hierarchicalDecomposition();
fprintf('Network decomposed into %d subsystems\n', length(sub_systems));

% ... continue with parameter estimation
```

## Main Components

### 1. FON_Graph_Class

The core class for fractional-order networks:

```matlab
% Create from state-space system
sys = ss(A, B, C, D);
fon = FON_Graph_Class(sys, alpha, G);

% Simulate the system
N = 100;  % Time steps
J = 50;   % Memory length
u = randn(m, N);  % Controllable inputs
w = randn(v, N);  % Disturbances
x0 = randn(n, 1); % Initial state
e = zeros(o, N);  % No measurement noise

fon.fsim(u, w, x0, e, J);

% Access simulation results
states = fon.x;   % State trajectories (n × N+1)
outputs = fon.y;  % Output trajectories (o × N)

% Perform decomposition
subsystems = fon.hierarchicalDecomposition();
```

### 2. Parameter Estimation Functions

**Zero-Input Estimation:**
```matlab
[est_alpha, est_A, est_error] = no_input_estimation_with_known_params(...
    x_data, J, A_known, A2est, alpha_known, alpha2est);
```

**Fractional Derivative Computation:**
```matlab
z = FON_Graph_Class.FON_z_sim(x, alpha, J);
```

### 3. Network Generation

```matlab
% Generate multiple random networks
n_trials = 50;
[A_trial, B_trial, C_trial, D_trial, G_trial, alpha_trial, x0_trial] = ...
    generateSparseNetwork(n, m, v, o, 'high', n_trials);
```

Sparsity levels: `'low'`, `'medium'`, `'high'`

## Key Parameters

### Simulation Parameters
- `N`: Number of time steps (recommended: 50-200)
- `J`: Memory length for fractional-order simulation (recommended: equal to N or N/2)
- `noise_sigma_w`: Process noise standard deviation (default: 0.5)
- `noise_sigma_e`: Measurement noise standard deviation (default: 0.25)

### Network Parameters
- `n`: Number of states (tested: 4-20)
- `m`: Number of controllable inputs (tested: 2-5)
- `v`: Number of uncontrollable inputs (tested: 1-3)
- `o`: Number of outputs (tested: 2-5)

## Computational Complexity

### Full System Estimation
- Time complexity: O(n³ × iterations)
- Memory: O(n² + n×N)

### Subsystem Decomposition
- Number of subsystems: typically ⌈n/3⌉ to ⌈n/2⌉
- Time complexity per subsystem: O(n_sub³ × iterations), where n_sub << n
- Overall speedup: 5-10× for networks with n > 15

## Performance Metrics

The code computes two main metrics:

1. **Mean Squared Error (MSE):**
   ```
   MSE = mean((y_true - y_est).^2)
   ```

2. **Pearson Correlation Coefficient:**
   ```
   ρ = corr(y_true', y_est')
   ```

Typical results:
- MSE: 0.01-0.10 (lower is better)
- ρ: 0.95-0.99 (higher is better)

## Troubleshooting

**Issue**: "Index exceeds array dimensions" during simulation  
**Solution**: Increase memory length `J` or reduce simulation time `N`

**Issue**: "Rank deficient" warning  
**Solution**: Increase the number of data points or check for structural identifiability issues

**Issue**: Poor estimation accuracy  
**Solution**:
- Increase SNR by reducing noise levels
- Increase simulation length `N`
- Check that network has sufficient input-output paths

**Issue**: FON_z_sim error about FON_Dtilde  
**Solution**: Make sure you're using the updated version with `FON_Graph_Class.FON_Dtilde`

## Citation

If you use this code in your research, please cite:

```bibtex
@article{yourname2025fon,
  title={Structural Identifiability in Fractional-Order Networks},
  author={Your Name and Coauthor Names},
  journal={Journal Name},
  year={2025},
  volume={XX},
  pages={XX--XX},
  doi={XX.XXXX/xxxxx}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions, issues, or collaboration:
- **Email:** [your.email@institution.edu]
- **GitHub Issues:** [Open an issue](https://github.com/yourusername/fractional-order-networks/issues)

## Acknowledgments

This work was supported by [funding agency/grant number]. We thank [collaborators/contributors] for their valuable feedback.

## Version History

- **v1.0** (January 2025): Initial release
  - Hierarchical decomposition algorithm
  - Two-stage parameter estimation
  - Comprehensive testing framework
