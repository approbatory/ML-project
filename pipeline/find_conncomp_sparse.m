function [cc, A] = find_conncomp_sparse(X, threshold, A2)    
% Finds connected components of columns of A wrt a corr. threshold
% A2 (optional) is an additional adjacency matrix as a multiplicative constraint
    
%     norms_X = sqrt(sum(X.^2,1));
%     X_norm = bsxfun(@rdivide,X,norms_X);
    %X_norm = zscore(X, 1, 1) / sqrt(size(X, 1));
    XCORR = spcorr(X);
    A = gather(XCORR > threshold);
    if exist('A2', 'var')
        A = A .* A2;
    end
    [s, memberships] = graphconncomp(sparse(A)); 
    cc = [];
    acc=0;
    for k = 1:s
        idx = find(memberships == k);
        if numel(idx) > 1
            % Find idx with most common neighbors
            [~, i] = max(sum(A(:, idx), 1));
            acc = acc + 1;
            % idx with most common neighbors goes last
            cc(acc).indices = [idx(setdiff(1:length(idx), i)), idx(i)];
        end
    end       
    
end
