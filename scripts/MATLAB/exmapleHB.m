addpath("src/MATLAB/");
addpath("src/MATLAB/utils/");
addpath("path_to_HierarchicalConsensus/");       % see README.md for details

% Example usage of the HCE method on Hierarchical Benchmark

% Generate a Hierarchical Benchmark model
p0 = 0.05;  % background probability of edges
p1 = 0.2;   % Probability of edges in the first level (second renormalized level)
p2 = 0.3;   % Probability of edges in the second level (first renormalized level)
p3 = 1 - p0 - p1 - p2; % Probability of edges in the fourth level (zeroth renormalized level)

[A, Sgtruth] = hierarchicalBenchmark(1000, [p0, p1, p2, p3]);
S = eventSamples(A, 100);
[Sc, Tree] = hierarchicalConsensus(S);

% Transform Tree to a linkage matrix (nTree)
nTree = buildLinkageFromTree(Sc, Tree);

if ~isempty(Tree)                           % Check if Tree is not empty (very rarely happens)
    % Zeroth renormalization level
    [labels, ~] = findHCELevel(nTree, [], 0);
    fprintf("Zeroth renormalization level: %f\n", AMI(labels, Sgtruth(:, end)));

    % First renormalization level
    [labels, ~] = findHCELevel(nTree, [], 1);
    fprintf("First renormalization level: %f\n", AMI(labels, Sgtruth(:, end - 1)));

    % Second renormalization level
    [labels, ~] = findHCELevel(nTree, [], 2);
    fprintf("Second renormalization level: %f\n", AMI(labels, Sgtruth(:, end - 2)));
end