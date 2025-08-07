classdef HNRG
    properties
        N
        R
        L
        Nh
        hierarchical_community_labels
        A
        plmax_f = false
    end
    
    methods
        function obj = HNRG(N, R, L, kav, rho, seed)
            if nargin < 4
                kav = 16;
            end
            if nargin < 5
                rho = 1;
            end
            if nargin == 6
                rng(seed);
            end
            
            obj.N = N;
            obj.R = R;
            obj.L = L;
            
            obj.Nh = N * (R + 1)^L;
            
            obj.hierarchical_community_labels = zeros(obj.Nh, 1);
            for l = 0:L-1
                nl = N * (R + 1)^(L - l - 1);
                communities = repmat((0:(R+1)^(l+1)-1)', 1, nl);
                communities = reshape(communities', [], 1);
                obj.hierarchical_community_labels = [obj.hierarchical_community_labels, communities];
            end

            % Sx computation
            Sx = zeros(1, L+1);
            Sx(1) = R * N * (R+1)^(L-1);
            for i = 1:L
                Sx(i+1) = N * (R+1)^(L - i);
            end
            
            % px computation
            px = containers.Map('KeyType', 'int64', 'ValueType', 'double');
            for i = 1:L
                px(i) = (rho^(L - i) / (1 + rho)^(L - i + 1)) * (kav / (Sx(i+1) - 1));
            end
            px(0) = (rho^L / (1 + rho)^L) * (kav / (Sx(1)));
            
            if px(L) > 1
                obj.plmax_f = true; 

            end
            
            % pij condensed pairwise edge probabilities
            condensed_N = obj.Nh * (obj.Nh - 1) / 2;
            pij = zeros(condensed_N, 1);
            
            x = 1;
            for i = 1:obj.Nh
                for j = (i+1):obj.Nh
                    ci = obj.hierarchical_community_labels(i, :);
                    cj = obj.hierarchical_community_labels(j, :);
                    match_level = find(ci == cj, 1, 'last') - 1; % adjust to 0-based level
                    pij(x) = px(match_level);
                    x = x + 1;
                end
            end
            
            % Generate adjacency matrix
            edges = pij > rand(size(pij));
            obj.A = squareform(edges); % MATLAB's squareform from Statistics and ML Toolbox
            obj.A = double(obj.A); % Convert logical to double
        end
    end
end