function gen = shufgen(X, ks, featmask)
if ~issparse(X)
    error('Can only shuffle sparse arrays');
end
ks = ks(:)';
K_vals = unique(ks);
[ii,jj,ss] = sp2cell(X);
[m,n] = size(X);
if length(featmask) ~= n
    error('featmask must have as many entries as features: has %d, needs %d', length(featmask), n);
end
gen = cell2sp(cellshuf(ii, ks, K_vals, featmask), jj, ss, m, n);
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
[i,j,s] = find(X);
Ns = diff([0; find(diff(j)); length(j)]);
ii = mat2cell(i, Ns);
jj = mat2cell(j, Ns);
ss = mat2cell(s, Ns);
end

function X = cell2sp(ii,jj,ss,m,n)
i = cell2mat(ii);
j = cell2mat(jj);
s = cell2mat(ss);
X = sparse(i,j,s,m,n);
end