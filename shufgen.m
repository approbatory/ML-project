function gen = shufgen(X, ks)
if ~issparse(X)
    error('Can only shuffle sparse arrays');
end
ks = ks(:)';
K_vals = unique(ks);
[ii,jj,ss] = sp2cell(X);
gen = @(featmask) cell2sp(cellshuf(ii, ks, K_vals, featmask), jj, ss);
end

function ii = cellshuf(ii, ks, K_vals, featmask)
for k = K_vals
    is_in_subset = ks == k;
    subset_inds = find(is_in_subset);
    num_in_class = sum(is_in_subset);
    for jj = 1:length(ii)
        if ~featmask(jj)
            continue;
        end
        i_mask = is_in_subset(ii{jj});
        num_to_shuf = sum(i_mask);
        ii{jj}(i_mask) = subset_inds(randperm(num_in_class, num_to_shuf));
    end
end
end

function [ii,jj,ss] = sp2cell(X)
[~,N] = size(X);
[i,j,s] = find(X);
Ns = zeros(N,1);
for j_ind = 1:N
    Ns(j_ind) = sum(j == j_ind);
end
ii = mat2cell(i, Ns);
jj = mat2cell(j, Ns);
ss = mat2cell(s, Ns);
end

function X = cell2sp(ii,jj,ss)
i = cell2mat(ii);
j = cell2mat(jj);
s = cell2mat(ss);
X = sparse(i,j,s);
end