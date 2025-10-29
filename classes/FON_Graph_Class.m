classdef FON_Graph_Class < handle
    properties
        % Graph structure, matrices, and parameters
        G               % Graph structure with X, U, W, Y, E fields

        % System matrices
        A               % State matrix
        B               % Input matrix
        G_matrix        % Uncontrollable input matrix
        C               % Output matrix
        D               % Feedthrough matrix (default: zero matrix)

        % System parameters
        alpha = []      % Fractional order vector
        n = 0           % Number of states
        m = 0           % Number of inputs
        v = 0           % Number of uncontrollable inputs
        o = 0           % Number of outputs

        % Decomposition results
        subsystems = {} % Cell array of subsystem graphs after decomposition

        % Simulation results
        x               % State trajectories matrix
        y               % Output trajectories matrix

        % Stability properties
        stable = false      % Stability flag
        eig_T = []         % System eigenvalues
    end


    methods
        function obj = FON_Graph_Class(varargin)
            % FON_Graph_Class Constructor
            % Usage options:
            %   1. FON_Graph_Class() - Creates a default model based on createCustomModel()
            %   2. FON_Graph_Class(G) - Initializes with custom graph G
            %   3. FON_Graph_Class(ss1, alpha) - Creates from state-space model and alpha values
            %
            % Inputs:
            %   G     - Custom graph structure (optional)
            %   ss1   - State space system (optional)
            %   alpha - Fractional order vector (optional)
            % Outputs:
            %   obj   - FON_Graph_Class object

            if nargin == 0
                % Create default model
                obj.G = obj.createCustomModel();

            elseif nargin == 1 && isstruct(varargin{1})
                % Initialize with custom graph
                obj.G = varargin{1};

                % Check if this is a subsystem structure (has X, U, W, Y, E fields)
                if isfield(obj.G, 'X') && isfield(obj.G, 'U') && ...
                        isfield(obj.G, 'W') && isfield(obj.G, 'Y') && isfield(obj.G, 'E')
                    % Extract dimensions from the graph
                    obj.n = length(obj.G.X);
                    obj.m = length(obj.G.U);
                    obj.v = length(obj.G.W);
                    obj.o = length(obj.G.Y);

                    % Extract system matrices and alpha from the graph structure
                    [obj.A, obj.B, obj.G_matrix, obj.C, obj.alpha] = obj.extractMatricesFromGraph(obj.G);

                    % Initialize feedthrough matrix to zeros
                    obj.D = zeros(obj.o, obj.m);
                else
                    error('Invalid graph structure provided');
                end

            elseif nargin >= 2 && isa(varargin{1}, 'ss')
                % Initialize from state-space model and alpha values
                ss1 = varargin{1};
                alpha = varargin{2};

                % Store system matrices
                obj.A = ss1.A;
                obj.B = ss1.B;
                obj.C = ss1.C;
                obj.D = ss1.D;
                obj.alpha = alpha;

                % Extract dimensions
                obj.n = size(obj.A, 1);
                obj.m = size(obj.B, 2);
                obj.o = size(obj.C, 1);

                % Assume no uncontrollable inputs if not specified
                if nargin >= 3
                    obj.G_matrix = varargin{3};
                    obj.v = size(obj.G_matrix, 2);
                else
                    obj.G_matrix = zeros(obj.n, 0);
                    obj.v = 0;
                end

                % Create graph representation
                obj.G = obj.createGraphFromMatrices(obj.A, obj.B, obj.G_matrix, obj.C, obj.alpha);
            end

            % % Only check stability if we have all necessary matrices and parameters
            % if ~isempty(obj.A) && ~isempty(obj.alpha) && length(obj.alpha) == size(obj.A, 1)
            %     try
            %         D_cells = obj.FON_Dtilde(zeros(1, obj.n), obj.alpha, 1);
            %         obj.eig_T = eig(obj.A - D_cells{1});
            %         obj.stable = all(abs(obj.eig_T) < 1);
            % 
            %         if obj.stable
            %             fprintf('The system is stable.\n');
            %         else
            %             fprintf('Warning: The system is unstable!\n');
            %         end
            %     catch e
            %         fprintf('Could not perform stability check: %s\n', e.message);
            %         obj.stable = false;
            %     end
            % end
        end

        function G_model = createCustomModel(obj)
            %CREATECUSTOMMODEL Creates a custom graph model
            %   This function creates a custom graph model for hierarchical decomposition.
            %   It defines a fractional-order network with specified parameters.

            %% Define system parameters
            n = 4;  % Number of states
            m = 2;  % Number of inputs
            v = 1;  % Number of uncontrollable inputs
            o = 2;  % Number of outputs

            % Define fractional-order exponents (α) - use distinct values
            %alpha = [0.3, 0.4, 0.5, 0.6];

            %% Define system matrices
            % A(θA) - state transition matrix
            % A = [
            %     0.3, 0.0, 0.0, 0.0;  % x1 has self-loop
            %     0.3, 0.4, 0.0, 0.0;  % x2 influenced by x1 and has self-loop
            %     0.0, 0.3, 0.5, 0.0;  % x3 influenced by x2 and has self-loop
            %     0.0, 0.0, 0.4, 0.6   % x4 influenced by x3 and has self-loop
            %     ];

            % New A matrix with more balanced dynamics
            A = [
                0.05, 0.00, 0.00, 0.00;  % Significantly reduced from 0.3
                0.10, 0.15, 0.00, 0.00;  % Reduced from 0.3, 0.4
                0.00, 0.10, 0.20, 0.00;  % Reduced from 0.3, 0.5
                0.00, 0.00, 0.10, 0.25   % Reduced from 0.4, 0.6
                ];

            % New alpha values (following pattern from your example)
            alpha = [0.15, 0.20, 0.25, 0.30];  %// Significantly reduced

            % B(θB) - input matrix
            B = zeros(n, m);
            B(1, 1) = 0.7;  % Input 1 affects state 1
            B(3, 2) = 0.8;  % Input 2 affects state 3

            % G(θG) - uncontrollable input matrix
            G_matrix = zeros(n, v);
            G_matrix(4, 1) = 0.5;  % Uncontrollable input affects state 4

            % C - output matrix
            C = zeros(o, n);
            C(1, 2) = 1;    % Output 1 measures state 2
            C(2, 4) = 1;    % Output 2 measures state 4

            %% Create the graph representation
            G_model = struct();

            % State vertices (numbered 1 to n)
            G_model.X = 1:n;

            % Input vertices (numbered from 101 to 100+m)
            G_model.U = 101:(100+m);

            % Uncontrollable input vertices (numbered from 201 to 200+v)
            G_model.W = 201:(200+v);

            % Output vertices (numbered from 301 to 300+o)
            G_model.Y = 301:(300+o);

            % Initialize edges
            edges = [];

            % Add edges from A matrix (state transitions) and fractional-order self-loops
            for i = 1:n
                for j = 1:n
                    if A(i, j) ~= 0 && i ~= j
                        % Edge from state j to state i (non-self-loop)
                        edges = [edges; j, i];
                    elseif (A(i, j) ~= 0 || alpha(i) ~= 0) && i == j
                        % Self-loop - only add once if either A(i,i) ≠ 0 or alpha(i) ≠ 0
                        edges = [edges; i, i];
                    end
                end
            end

            % Add edges from B matrix (input to state)
            for i = 1:n
                for j = 1:m
                    if B(i, j) ~= 0
                        edges = [edges; 100+j, i];  % Edge from input j to state i
                    end
                end
            end

            % Add edges from G matrix (uncontrollable input to state)
            for i = 1:n
                for j = 1:v
                    if G_matrix(i, j) ~= 0
                        edges = [edges; 200+j, i];  % Edge from uncontrollable input j to state i
                    end
                end
            end

            % Add edges from C matrix (state to output)
            for i = 1:o
                for j = 1:n
                    if C(i, j) ~= 0
                        edges = [edges; j, 300+i];  % Edge from state j to output i
                    end
                end
            end

            % Assign edges to the graph
            G_model.E = edges;

            % Store these matrices in the object for future use
            obj.A = A;
            obj.B = B;
            obj.G_matrix = G_matrix;
            obj.C = C;
            obj.alpha = alpha;
            obj.n = n;
            obj.m = m;
            obj.v = v;
            obj.o = o;
        end

        function G_model = createGraphFromMatrices(obj, A, B, G_matrix, C, alpha)
            %CREATEGRAPHFROMMATRICES Creates a graph representation from system matrices
            %   This function creates a graph model based on provided system matrices

            % Get dimensions
            n = size(A, 1);
            m = size(B, 2);
            v = size(G_matrix, 2);
            o = size(C, 1);

            % Create graph structure
            G_model = struct();

            % State vertices (numbered 1 to n)
            G_model.X = 1:n;

            % Input vertices (numbered from 101 to 100+m)
            G_model.U = 101:(100+m);

            % Uncontrollable input vertices (numbered from 201 to 200+v)
            G_model.W = 201:(200+v);

            % Output vertices (numbered from 301 to 300+o)
            G_model.Y = 301:(300+o);

            % Initialize edges
            edges = [];

            % Add edges from A matrix (state transitions) and fractional-order self-loops
            for i = 1:n
                for j = 1:n
                    if A(i, j) ~= 0 && i ~= j
                        % Edge from state j to state i (non-self-loop)
                        edges = [edges; j, i];
                    elseif (A(i, j) ~= 0 || alpha(i) ~= 0) && i == j
                        % Self-loop - only add once if either A(i,i) ≠ 0 or alpha(i) ≠ 0
                        edges = [edges; i, i];
                    end
                end
            end

            % Add edges from B matrix (input to state)
            for i = 1:n
                for j = 1:m
                    if B(i, j) ~= 0
                        edges = [edges; 100+j, i];  % Edge from input j to state i
                    end
                end
            end

            % Add edges from G matrix (uncontrollable input to state)
            for i = 1:n
                for j = 1:v
                    if G_matrix(i, j) ~= 0
                        edges = [edges; 200+j, i];  % Edge from uncontrollable input j to state i
                    end
                end
            end

            % Add edges from C matrix (state to output)
            for i = 1:o
                for j = 1:n
                    if C(i, j) ~= 0
                        edges = [edges; j, 300+i];  % Edge from state j to output i
                    end
                end
            end

            % Assign edges to the graph
            G_model.E = edges;
        end

        function [A, B, G_matrix, C, alpha] = extractMatricesFromGraph(obj, G)
            %EXTRACTMATRICESFROMGRAPH Extracts system matrices from a graph structure
            %   This function reconstructs the system matrices from a graph representation

            % Get dimensions
            n = length(G.X);
            m = length(G.U);
            v = length(G.W);
            o = length(G.Y);

            % Initialize matrices
            A = zeros(n, n);
            B = zeros(n, m);
            G_matrix = zeros(n, v);
            C = zeros(o, n);
            alpha = zeros(1, n);  % Initialize alpha vector

            % Set default alpha values if needed
            for i = 1:n
                % Check if state i has a self-loop in the graph
                if any(all(G.E == [i, i], 2))
                    alpha(i) = 0.7; % Default value, can be adjusted
                end
            end

            % Process edges to reconstruct matrices
            for i = 1:size(G.E, 1)
                from_node = G.E(i, 1);
                to_node = G.E(i, 2);

                % State to state edges (A matrix)
                if ismember(from_node, G.X) && ismember(to_node, G.X)
                    A(to_node, from_node) = 0.5; % Default weight if not specified
                end

                % Input to state edges (B matrix)
                if ismember(from_node, G.U) && ismember(to_node, G.X)
                    input_idx = from_node - 100;
                    B(to_node, input_idx) = 0.5; % Default weight if not specified
                end

                % Uncontrollable input to state edges (G matrix)
                if ismember(from_node, G.W) && ismember(to_node, G.X)
                    w_idx = from_node - 200;
                    G_matrix(to_node, w_idx) = 0.5; % Default weight if not specified
                end

                % State to output edges (C matrix)
                if ismember(from_node, G.X) && ismember(to_node, G.Y)
                    output_idx = to_node - 300;
                    C(output_idx, from_node) = 1; % Default weight
                end
            end
        end

        function subsystem_objects = hierarchicalDecomposition(obj)
            %HIERARCHICALDECOMPOSITION Decomposes a system graph into subsystems
            %   This function implements a hierarchical decomposition algorithm to divide
            %   a system graph into subsystems, analyzing one output at a time and finding
            %   paths connecting each output to an input with minimal overlap.
            %
            %   Output:
            %   subsystem_objects - Cell array of FON_Graph_Class objects, each representing a subsystem

            % First, perform the original decomposition to get subsystem structures
            subsystems = obj.decomposeGraphIntoSubsystems(obj.G);

            % Store the subsystems in the object
            obj.subsystems = subsystems;

            % Create FON_Graph_Class objects for each subsystem
            subsystem_objects = cell(1, length(subsystems));

            for i = 1:length(subsystems)
                % Create a new FON_Graph_Class object initialized with the subsystem graph
                subsystem_objects{i} = FON_Graph_Class(subsystems{i});
            end
        end


        function subsystems = decomposeGraphIntoSubsystems(obj, G)
            % DECOMPOSEGRAPHINTOSUBSYSTEMS Backward decomposition following Algorithm 1
            % Starts from each output and traces backwards to find all reachable nodes
            %
            % Input:
            %   G - Graph struct with fields: X, U, W, Y, E
            %
            % Output:
            %   subsystems - Cell array of subgraph structures

            subsystems = {};

            % Get all output vertices and sort them
            output_vertices = sort(G.Y);

            % Process each output in order (y1, y2, y3, ...)
            for i = 1:length(output_vertices)
                out_vertex = output_vertices(i);

                % Perform backward reachability from this output
                reachable_vertices = performBackwardReach(G, out_vertex);

                % Create subgraph with reachable vertices
                subG = createSubgraph(G, reachable_vertices);

                % Store subsystem
                subsystems{end+1} = subG;
            end

            % Nested helper functions
            function reachable_vertices = performBackwardReach(graph, start_vertex)
                % Perform backward reachability analysis from start_vertex
                % Includes bidirectional connections between states

                reachable_vertices = start_vertex;
                vertices_to_explore = start_vertex;
                explored_vertices = [];

                while ~isempty(vertices_to_explore)
                    current_vertex = vertices_to_explore(1);
                    vertices_to_explore(1) = [];

                    if ismember(current_vertex, explored_vertices)
                        continue;
                    end
                    explored_vertices = [explored_vertices, current_vertex];

                    % Find predecessors (vertices with edges TO current_vertex)
                    predecessors = graph.E(graph.E(:,2) == current_vertex, 1);

                    % Find successors (vertices with edges FROM current_vertex)
                    successors = graph.E(graph.E(:,1) == current_vertex, 2);

                    % For state nodes, include bidirectional connections
                    if ismember(current_vertex, graph.X)
                        % Check successors for bidirectional edges
                        for succ = successors'
                            if ismember(succ, graph.X) && ~ismember(succ, reachable_vertices)
                                % Check if there's an edge back (bidirectional)
                                if any(graph.E(:,1) == succ & graph.E(:,2) == current_vertex)
                                    reachable_vertices = [reachable_vertices, succ];
                                    vertices_to_explore = [vertices_to_explore, succ];
                                end
                            end
                        end
                    end

                    % Add all predecessors
                    for pred = predecessors'
                        if ~ismember(pred, reachable_vertices)
                            reachable_vertices = [reachable_vertices, pred];
                            vertices_to_explore = [vertices_to_explore, pred];
                        end
                    end
                end

                % Sort vertices
                reachable_vertices = sort(reachable_vertices);
            end

            function subG = createSubgraph(graph, vertices)
                % Extract subgraph containing only the specified vertices

                subG = struct();

                % Separate vertices by type
                subG.X = vertices(ismember(vertices, graph.X));
                subG.U = vertices(ismember(vertices, graph.U));
                subG.W = vertices(ismember(vertices, graph.W));
                subG.Y = vertices(ismember(vertices, graph.Y));

                % Extract edges that connect vertices in the subgraph
                subG.E = [];
                for i = 1:size(graph.E, 1)
                    from_v = graph.E(i, 1);
                    to_v = graph.E(i, 2);

                    % Include edge if both vertices are in the subgraph
                    if ismember(from_v, vertices) && ismember(to_v, vertices)
                        subG.E = [subG.E; from_v, to_v];
                    end
                end
            end
        end


        % % no overlapping decomposition
        % function subsystems = decomposeGraphIntoSubsystems(obj, G)
        %     % DECOMPOSEGRAPHINTOSUBSYSTEMS Backward decomposition following Algorithm 1
        %     % Starts from each output and traces backwards to find all reachable nodes
        %     % States already claimed by previous subsystems act as "cut points"
        %     %
        %     % Input:
        %     %   G - Graph struct with fields: X, U, W, Y, E
        %     %
        %     % Output:
        %     %   subsystems - Cell array of subgraph structures
        % 
        %     subsystems = {};
        %     claimed_states = [];  % Track states already assigned to previous subsystems
        % 
        %     % Get all output vertices and sort them
        %     output_vertices = sort(G.Y);
        % 
        %     % Process each output in order (y1, y2, y3, ...)
        %     for i = 1:length(output_vertices)
        %         out_vertex = output_vertices(i);
        % 
        %         % Perform backward reachability from this output, stopping at claimed states
        %         reachable_vertices = performBackwardReach(G, out_vertex, claimed_states);
        % 
        %         % Create subgraph with reachable vertices
        %         subG = createSubgraph(G, reachable_vertices);
        % 
        %         % Store subsystem
        %         subsystems{end+1} = subG;
        % 
        %         % Update claimed states with new states from this subsystem
        %         new_states = reachable_vertices(ismember(reachable_vertices, G.X));
        %         claimed_states = [claimed_states, new_states];
        %     end
        % 
        %     % Nested helper functions
        %     function reachable_vertices = performBackwardReach(graph, start_vertex, claimed_states)
        %         % Perform backward reachability analysis from start_vertex
        %         % STOP at any state that has already been claimed by a previous subsystem
        %         % Includes bidirectional connections between states
        % 
        %         reachable_vertices = start_vertex;
        %         vertices_to_explore = start_vertex;
        %         explored_vertices = [];
        % 
        %         while ~isempty(vertices_to_explore)
        %             current_vertex = vertices_to_explore(1);
        %             vertices_to_explore(1) = [];
        % 
        %             if ismember(current_vertex, explored_vertices)
        %                 continue;
        %             end
        %             explored_vertices = [explored_vertices, current_vertex];
        % 
        %             % If current vertex is a state that's already claimed, STOP here (cut point)
        %             if ismember(current_vertex, graph.X) && ismember(current_vertex, claimed_states)
        %                 continue;  % Don't explore beyond this claimed state
        %             end
        % 
        %             % Find predecessors (vertices with edges TO current_vertex)
        %             predecessors = graph.E(graph.E(:,2) == current_vertex, 1);
        % 
        %             % Find successors (vertices with edges FROM current_vertex)
        %             successors = graph.E(graph.E(:,1) == current_vertex, 2);
        % 
        %             % For state nodes, include bidirectional connections
        %             if ismember(current_vertex, graph.X)
        %                 % Check successors for bidirectional edges
        %                 for succ = successors'
        %                     if ismember(succ, graph.X) && ~ismember(succ, reachable_vertices)
        %                         % Don't add if it's already claimed
        %                         if ismember(succ, claimed_states)
        %                             continue;
        %                         end
        % 
        %                         % Check if there's an edge back (bidirectional)
        %                         if any(graph.E(:,1) == succ & graph.E(:,2) == current_vertex)
        %                             reachable_vertices = [reachable_vertices, succ];
        %                             vertices_to_explore = [vertices_to_explore, succ];
        %                         end
        %                     end
        %                 end
        %             end
        % 
        %             % Add all predecessors (but skip claimed states)
        %             for pred = predecessors'
        %                 if ~ismember(pred, reachable_vertices)
        %                     % If predecessor is a claimed state, don't add it
        %                     if ismember(pred, graph.X) && ismember(pred, claimed_states)
        %                         continue;
        %                     end
        % 
        %                     reachable_vertices = [reachable_vertices, pred];
        %                     vertices_to_explore = [vertices_to_explore, pred];
        %                 end
        %             end
        %         end
        % 
        %         % Sort vertices
        %         reachable_vertices = sort(reachable_vertices);
        %     end
        % 
        %     function subG = createSubgraph(graph, vertices)
        %         % Extract subgraph containing only the specified vertices
        % 
        %         subG = struct();
        % 
        %         % Separate vertices by type
        %         subG.X = vertices(ismember(vertices, graph.X));
        %         subG.U = vertices(ismember(vertices, graph.U));
        %         subG.W = vertices(ismember(vertices, graph.W));
        %         subG.Y = vertices(ismember(vertices, graph.Y));
        % 
        %         % Extract edges that connect vertices in the subgraph
        %         subG.E = [];
        %         for i = 1:size(graph.E, 1)
        %             from_v = graph.E(i, 1);
        %             to_v = graph.E(i, 2);
        % 
        %             % Include edge if both vertices are in the subgraph
        %             if ismember(from_v, vertices) && ismember(to_v, vertices)
        %                 subG.E = [subG.E; from_v, to_v];
        %             end
        %         end
        %     end
        % end
        % 

        function adj_list = createAdjacencyList(obj, G)
            % Helper function to create an adjacency list for the graph
            % Find maximum vertex ID to determine size
            max_id = max([max(G.X), max(G.U), max(G.W), max(G.Y)]);
            adj_list = cell(max_id, 1);

            % Add all edges to the adjacency list
            for i = 1:size(G.E, 1)
                from = G.E(i, 1);
                to = G.E(i, 2);

                % Initialize if needed
                if isempty(adj_list{from})
                    adj_list{from} = [];
                end

                % Add neighbor
                adj_list{from} = [adj_list{from}, to];
            end
        end


        function paths = findPathsBFS(obj, adj_list, start, target, all_states, covered_states)
            % BFS-based path finding with strong preference for non-overlapping paths
            paths = {};

            % Create a queue for BFS
            queue = cell(1, 1000);  % Pre-allocate for efficiency
            queue{1} = [start];
            queue_start = 1;
            queue_end = 1;

            % Track visited nodes to avoid cycles
            max_node = length(adj_list);
            visited = false(1, max_node);

            % BFS with path tracking
            while queue_start <= queue_end
                % Get next path from queue
                current_path = queue{queue_start};
                queue_start = queue_start + 1;

                % Get the last node in the path
                current = current_path(end);

                % If we've reached the target, add this path to our results
                if current == target
                    paths{end+1} = current_path;

                    % Limit the number of paths we find
                    if length(paths) >= 10
                        break;
                    end
                    continue;
                end

                % Mark as visited
                visited(current) = true;

                % Get neighbors
                if current <= length(adj_list) && ~isempty(adj_list{current})
                    neighbors = adj_list{current};

                    % First prioritize neighbors that aren't covered yet
                    if ~isempty(neighbors)
                        state_neighbors = neighbors(ismember(neighbors, all_states));

                        % Sort neighbors: uncovered first, then covered
                        uncovered = state_neighbors(~ismember(state_neighbors, covered_states));
                        covered = state_neighbors(ismember(state_neighbors, covered_states));
                        sorted_neighbors = [uncovered, covered, setdiff(neighbors, state_neighbors)];

                        % Process each neighbor
                        for next = sorted_neighbors
                            % Skip if already visited in this path
                            if any(current_path == next) || visited(next)
                                continue;
                            end

                            % Create new path and add to queue
                            new_path = [current_path, next];
                            queue_end = queue_end + 1;
                            queue{queue_end} = new_path;

                            % Limit path length to prevent excessive memory use
                            if length(new_path) > 10
                                break;
                            end
                        end
                    end
                end
            end
        end


        function visualizeGraph(obj)
            % VISUALIZEGRAPH Visualizes the system graph with larger elements

            figure('Position', [100, 100, 1000, 800]);  % Larger figure window

            G_vis = digraph();

            % Add nodes with proper labels
            for x = obj.G.X
                G_vis = addnode(G_vis, ['x' num2str(x)]);
            end
            for u = obj.G.U
                G_vis = addnode(G_vis, ['u' num2str(u-100)]);
            end
            for w = obj.G.W
                G_vis = addnode(G_vis, ['w' num2str(w-200)]);
            end
            for y = obj.G.Y
                G_vis = addnode(G_vis, ['y' num2str(y-300)]);
            end

            % Add edges
            for i = 1:size(obj.G.E, 1)
                from_v = obj.G.E(i, 1);
                to_v = obj.G.E(i, 2);

                % Convert to proper node name
                if ismember(from_v, obj.G.X)
                    from_label = ['x' num2str(from_v)];
                elseif ismember(from_v, obj.G.U)
                    from_label = ['u' num2str(from_v-100)];
                elseif ismember(from_v, obj.G.W)
                    from_label = ['w' num2str(from_v-200)];
                else
                    from_label = ['y' num2str(from_v-300)];
                end

                if ismember(to_v, obj.G.X)
                    to_label = ['x' num2str(to_v)];
                elseif ismember(to_v, obj.G.U)
                    to_label = ['u' num2str(to_v-100)];
                elseif ismember(to_v, obj.G.W)
                    to_label = ['w' num2str(to_v-200)];
                else
                    to_label = ['y' num2str(to_v-300)];
                end

                G_vis = addedge(G_vis, from_label, to_label);
            end

            % Plot with increased size
            h = plot(G_vis, 'Layout', 'layered');
            title('Fractional-Order Network Graph', 'FontSize', 24, 'FontWeight', 'bold');

            % Increase line width (edges)
            h.LineWidth = 3.5;

            % Increase marker size (nodes)
            h.MarkerSize = 20;

            % Increase arrow size
            h.ArrowSize = 15;

            % Store original node labels
            node_labels = h.NodeLabel;

            % Hide the default node labels
            h.NodeLabel = {};

            % Get node positions
            x_pos = h.XData;
            y_pos = h.YData;

            % Add custom text labels with offset
            offset = 0.15;  % Adjust this value to move labels further/closer
            for i = 1:length(node_labels)
                text(x_pos(i), y_pos(i) + offset, node_labels{i}, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 18, ...
                    'FontWeight', 'bold', ...
                    'BackgroundColor', 'white', ...
                    'EdgeColor', 'none', ...
                    'Margin', 2);
            end

            % Assign node colors correctly
            state_nodes = findnode(G_vis, cellfun(@(x) ['x' num2str(x)], num2cell(obj.G.X), 'UniformOutput', false));
            input_nodes = findnode(G_vis, cellfun(@(x) ['u' num2str(x-100)], num2cell(obj.G.U), 'UniformOutput', false));
            uncont_nodes = findnode(G_vis, cellfun(@(x) ['w' num2str(x-200)], num2cell(obj.G.W), 'UniformOutput', false));
            output_nodes = findnode(G_vis, cellfun(@(x) ['y' num2str(x-300)], num2cell(obj.G.Y), 'UniformOutput', false));

            highlight(h, state_nodes, 'NodeColor', 'b');
            highlight(h, input_nodes, 'NodeColor', 'g');
            highlight(h, uncont_nodes, 'NodeColor', 'r');
            highlight(h, output_nodes, 'NodeColor', 'm');

            % Create a legend with larger font
            leg = legend('States', 'Inputs', 'Uncontrollable Inputs', 'Outputs', 'Location', 'best');
            leg.FontSize = 24;
        end


        function fsim(obj, u, w, x0, e, J)
            % Simulate the FON system dynamics
            % Inputs:
            %   u  - Input sequence
            %   w  - Uncontrollable input sequence
            %   x0 - Initial state
            %   e  - Additive noise in output
            %   J  - Memory length parameter
            % Outputs:
            %   Updates obj.x - State trajectories matrix (n × N)
            %   Updates obj.y - Output trajectories matrix (o × N)

            validateattributes(u, {'numeric'}, {'2d'}, 'fsim', 'u');
            validateattributes(w, {'numeric'}, {'2d'}, 'fsim', 'w');
            validateattributes(x0, {'numeric'}, {'vector', 'numel', obj.n}, 'fsim', 'x0');
            validateattributes(e, {'numeric'}, {'2d'}, 'fsim', 'e');
            validateattributes(J, {'numeric'}, {'positive', 'integer', 'scalar'}, 'fsim', 'J');

            % Initialize simulation arrays
            N = size(u, 2);
            obj.x = zeros(obj.n, N+1);
            obj.y = zeros(obj.o, N);

            % Set initial state
            obj.x(:,1) = x0;

            % Pre-compute D matrices for the entire simulation
            D_cells = obj.FON_Dtilde(zeros(size(obj.alpha)), obj.alpha, min(N,J));  % CHANGED: Added min(N,J)

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

                % State update
                obj.x(:,i) = obj.A * obj.x(:,i-1) + obj.B * u(:,i-1) + obj.G_matrix * w(:,i-1) ...
                    - frac_sum;

                % Output equation with D matrix (direct feedthrough)
                obj.y(:,i-1) = obj.C * obj.x(:,i) + obj.D * u(:,i-1) + e(:,i-1);  % CHANGED: Added obj.D * u(:,i-1)
            end
        end

    end

    methods (Static)
        function psi = FON_psi(alpha, j)
            % Compute psi values for a scalar alpha
            % Inputs:
            %   alpha - Scalar fractional order
            %   j     - Number of iterations
            % Outputs:
            %   psi   - Column vector ((j+1) × 1) containing psi values
            %          starting from psi(alpha,0)

            validateattributes(alpha, {'numeric'}, {'scalar'});

            psi(1,1) = 1;
            for i = 1:j
                psi(i+1,1) = psi(i,1) * ((i-1-alpha)/i);
            end
        end

        function D = FON_Dtilde(alpha0, delta_alpha, k)
            % Generate Dk matrices using the simplified formulation
            % Inputs:
            %   alpha0      - Vector of initial alpha values [alpha0_1, alpha0_2, ...]
            %   delta_alpha - Vector of alpha changes [delta_alpha_1, delta_alpha_2, ...]
            %   k           - Maximum time index
            % Outputs:
            %   D           - Cell array {k} where D{j} is the diagonal matrix D at time j

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
                    psi_alpha_plus_delta(:,i) = FON_Graph_Class.FON_psi(delta_alpha(i), k);
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
                    psi_alpha(:,i) = FON_Graph_Class.FON_psi(alpha0(i), k);
                    psi_alpha_plus_delta(:,i) = FON_Graph_Class.FON_psi(alpha0(i) + delta_alpha(i), k);
                end

                % For each time step j, compute D{j} = Dtilde(0,alpha+delta_alpha,j) - Dtilde(0,alpha,j)
                for j = 1:k
                    D{j} = diag(psi_alpha_plus_delta(j+1,:) - psi_alpha(j+1,:));
                end
            end
        end
    end



end