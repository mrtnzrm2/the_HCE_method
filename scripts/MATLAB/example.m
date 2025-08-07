addpath("src/MATLAB/");
addpath("src/MATLAB/utils/");

% Example usage of the HCE method on a HNRG model
% Generate a HNRG model

N = 10 ; % number of nodes
R = 3 ; % Branching factor
L = 3 ; % Number of levels
kav = 16 ; % Average degree
rho = 1 ; % Cohesivness

% Generate the HNRG model
G = HNRG(N, R, L, kav, rho, "shuffle");
D = compute_dissimilarity_matrix(G.A);
H = linkage(squareform(D), 'average');

% Zeroth renormalization level
[labels, ~] = findHCELevel(H, [], 0);

fprintf("Zeroth renormalization level: %f\n", AMI(labels, G.hierarchical_community_labels(:, end) + 1));

% First renormalization level
[labels, ~] = findHCELevel(H, [], 1);
fprintf("First renormalization level: %f\n", AMI(labels, G.hierarchical_community_labels(:, end-1) + 1));

% Second renormalization level
[labels, ~] = findHCELevel(H, [], 2);
fprintf("Second renormalization level: %f\n", AMI(labels, G.hierarchical_community_labels(:, end-2) + 1));
