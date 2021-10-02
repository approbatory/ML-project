function X_shuf = fshuffle(X, ks, auto)
if ~exist('auto', 'var')
    auto = false;
end

X_shuf = nan(size(X));
K_vals = unique(ks);
K = numel(K_vals);

P = size(X, 2);

for k_ix = 1:K
    K_val = K_vals(k_ix);
    [pass_s, pass_e] = get_passes(ks == K_val);
    n_passes = numel(pass_s);
    for p_ix = 1:P
        for pass_ix = 1:n_passes
            if auto
                other_pass = pass_ix;
            else
                other_pass = randi(n_passes);
            end
            other_pass_data = X(pass_s(other_pass):pass_e(other_pass),p_ix);
            current_pass_length = pass_e(pass_ix) - pass_s(pass_ix) + 1;
            resampled_pass = take_samples(current_pass_length, other_pass_data);
            X_shuf(pass_s(pass_ix):pass_e(pass_ix), p_ix) = resampled_pass;
        end
    end
end

assert(~any(isnan(X_shuf(:))));
end

function [pass_begin_idx, pass_end_idx] = get_passes(b)
b = b(:);
b_aug = [0; b; 0];
diff_b_aug = diff(b_aug);
pass_begin_idx = find(diff_b_aug == 1);
pass_end_idx = find(diff_b_aug == -1) - 1;
assert(numel(pass_begin_idx) == numel(pass_end_idx));
end

function s = take_samples(n, given)
n_given = numel(given);
s_idx = randi(n_given, n, 1);
s = given(s_idx);
end