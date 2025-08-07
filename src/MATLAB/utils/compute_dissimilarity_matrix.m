function DistA = compute_dissimilarity_matrix(A, f)
% Computes the dissimilarity matrix for a given adjacency matrix A.
%
% Parameters
% ----------
% A : matrix
%     Adjacency matrix of the graph.
% f : function handle, optional
%     Function to compute the dissimilarity between two nodes. Default is
%     @graph_cosine_similarity.
%
% Returns
% -------
% DistA : vector
%     Dissimilarity matrix, a condensed 1D array of distances between nodes.

    if nargin < 2
        f = @graph_cosine_similarity;
    end

    if ~ismatrix(A) || size(A,1) ~= size(A,2)
        error('Input must be a square adjacency matrix.');
    end

    nodes = size(A, 1);
    DistA = zeros(nodes*(nodes-1)/2, 1);
    e = 1;
    for i = 1:nodes
        for j = i+1:nodes
            DistA(e) = f(A(i,:), A(j,:), i, j);
            e = e + 1;
        end
    end

    DistA = min(max(DistA, -1), 1);    % Equivalent to np.clip
    DistA = sqrt(2 * (1 - DistA));
    DistA(isnan(DistA)) = 2;
end