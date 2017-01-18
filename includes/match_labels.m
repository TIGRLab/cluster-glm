function [mappings] = match_labels(a, b);
    % match_labels(a, b)
    %
    % returns mappings from a to b that minimizes the distance between the 
    % input label vectors.
    % inputs must have the same number of labels (k) and same length.
    %
    % For details see section 2.2 of:
    %    Lange T, Roth V, Braun M, Buhmann J. 2004. Stability-based validation
    %    of clustering solutions. Neural computation 16, 1299-1323.  
    %   
    % a = [1,1,2,3,3,4,4,4,2];
    % b = [2,2,3,1,1,4,4,4,3];
    % optimal: 1 => 2, 2 => 3, 3 => 1, 4 => 4
    %
    % mappings =
    %      1     2
    %      2     3
    %      3     1
    %      4     4
    %
    % Inspired by http://things-about-r.tumblr.com/post/36087795708/matching-clustering-solutions-using-the-hungarian
    
    ids_a = unique(a);
    ids_b = unique(b);
    n = length(ids_a);
    matrix = zeros(n); % distance matrix (between labels)

    % check the inputs
    if length(ids_a) ~= length(ids_b) || length(a) ~= length(b);
        disp('The length of labels a + b should be the same, and')
        disp('they should both consist of k labels.')
        
        return
    end

    % construct the distance matrix
    for x = 1:n;
        for y = 1:n;
            idx_a = find(a == x);
            idx_b = find(b == y);
            n_int = length(intersect(idx_a, idx_b));
            % distance = sum(# in cluster) - 2*sum(# in intersection)
            matrix(x,y) = (length(idx_a) + length(idx_b) - 2*n_int);
        end
    end

    % permute labels using minimum weighted bipartite matching
    matrix = hungarian(matrix);

    % get mappings
    mappings = zeros(n,1);
    for x = 1:n;
        mappings(x) = find(matrix(x, :) == 1);
    end
end
