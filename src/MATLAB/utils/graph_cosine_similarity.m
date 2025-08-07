function sim = graph_cosine_similarity(u, v, i, j)
% Computes the cosine similarity between two vectors u and v, excluding the
% i-th and j-th elements respectively.
%
% Parameters
% ----------
% u : vector
%     First vector.
% v : vector
%     Second vector.
% i : int
%     Index of the node i in the vector u.
% j : int
%     Index of the node j in the vector v.
%
% Returns
% -------
% sim : float
%     Cosine similarity between u and v, excluding the specified elements.

    if ~isvector(u) || ~isvector(v)
        error('Both u and v must be 1D arrays.');
    end

    uc = u;
    vc = v;

    % MATLAB uses 1-based indexing, so ensure i and j are correct
    uc(j) = 0;
    vc(i) = 0;

    uj = u(j);
    vi = v(i);

    uv = dot(uc, vc) + uj*vi + u(i)*v(j);
    uu = dot(u, u);
    vv = dot(v, v);

    if uu == 0 || vv == 0
        sim = -1;
    else
        sim = uv / sqrt(uu * vv);
    end
end