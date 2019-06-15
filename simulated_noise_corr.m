%simulated shuf, unshuf, diag
n_series = 1:1:500;
unshuf_dprime2 = zeros(size(n_series));
shuf_dprime2 = zeros(size(n_series));
diag_dprime2 = zeros(size(n_series));
for n_i = 1:numel(n_series)
    n = n_series(n_i);
    m_vec_ = rand(n,1)*2-1;%rand(n,1);
    Sigma_mat_ = eye(n);
    for i = 1:n
        for j = 1:i-1
            Sigma_mat_(i,j) = 0.1*rand;
            Sigma_mat_(j,i) = Sigma_mat_(i,j);
        end
    end
    D_mat_ = diag(diag(Sigma_mat_));
    
    %diag_result_  = (m_vec_' * (D_mat_ \ Sigma_mat_) * (D_mat_ \ m_vec_)) / (m_vec_' * (D_mat_ \ (D_mat_ \ m_vec_)));
    %unshuf_result_ = (m_vec_' * (Sigma_mat_ \ m_vec_)) / (m_vec_' * (Sigma_mat_ \ (Sigma_mat_ \ m_vec_)));
    %shuf_result_ = (m_vec_' * (D_mat_ \ m_vec_)) / (m_vec_' * (D_mat_ \ (D_mat_ \ m_vec_)));
    
    wopt = Sigma_mat_ \ m_vec_;
    unshuf_dprime2(n_i) = m_vec_' * wopt;
    
    wdiag = D_mat_ \ m_vec_;
    shuf_dprime2(n_i) = m_vec_' * wdiag;
    
    diag_dprime2(n_i) = unshuf_dprime2(n_i).^2 / (wdiag' * Sigma_mat_ * wdiag);
    
    
    %diag_result_(n_i)  = m_vec_' * inv(D_mat_) * Sigma_mat_ * inv(D_mat_) * m_vec_ / (m_vec_' * inv(D_mat_) * inv(D_mat_) * m_vec_);
    %unshuf_result_(n_i) = m_vec_' * inv(Sigma_mat_) * m_vec_ / (m_vec_' * inv(Sigma_mat_) * inv(Sigma_mat_) * m_vec_);
    %shuf_result_(n_i) = m_vec_' * inv(D_mat_) * m_vec_ / (m_vec_' * inv(D_mat_) * inv(D_mat_) * m_vec_);
    
    %unshuf_dprime2(n_i) = m_vec_' * inv(Sigma_mat_) * m_vec_;
    %shuf_dprime2(n_i) = m_vec_' * inv(D_mat_) * m_vec_;
    %diag_dprime2(n_i) = unshuf_dprime2(n_i).^2 / (m_vec_' * inv(D_mat_) * Sigma_mat_ * inv(D_mat_) * m_vec_);
    
    %fprintf('n=%d, Diag noise: %.2e, Unshuf noise: %.2e, Shuf noise: %.2e\n', n, diag_result_(n_i), unshuf_result_(n_i), shuf_result_(n_i));
    fprintf('n=%d, Diag dp2: %.2e, Unshuf dp2: %.2e, Shuf dp2: %.2e\n', n, diag_dprime2(n_i), unshuf_dprime2(n_i), shuf_dprime2(n_i));
end
%%
figure;
plot(n_series, unshuf_dprime2, 'b');
hold on;
plot(n_series, shuf_dprime2, 'r');
plot(n_series, diag_dprime2, 'm');