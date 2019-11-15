function [ Xdpca, Xphi, margNums ] = anova_decomp( X, meta )
%ANOVA_DECOMP Analysis-of-variance decomposition of data tensor
%   Xphi = ANOVA_DECOMP(X, meta)

%% make Xdpca -> (neurons x time x start x day) array of averages

% data dimensions
[N,T,K] = size(X);

% mean center
X = X - repmat(mean(X(:,:),2), 1, T, K);

% get indices for starts
s = zeros(K,1);
s(strcmp(meta.start,'east')) = 1;
s(strcmp(meta.start,'west')) = 2;

% get indices for day
days = unique(meta.day);
ndays = length(days);
d = zeros(K,1);
for d_ = 1:ndays
    d(meta.day == days(d_)) = d_;
end

% find averages
Xdpca = nan(N,T,2,ndays);
for s_ = 1:2
    for d_ = 1:ndays
        idx = (s==s_) & (d==d_);
        Xdpca(:,:,s_,d_) = mean(X(:,:,idx),3);
    end
end

%% marginalize
combinedParams = {{1}, {2, [1 2]}, {3, [1 3]}, {[2 3], [1 2 3]}};
[Xphi, margNums] = dpca_marginalize(Xdpca,'combinedParams',combinedParams);


%% visualize
neurons = randperm(N);
for n = neurons(1:10)
    figure()
    subplot(4,1,1); hold on
    plot(Xphi{1}(n,:,1,1))
    title('trial average')
    
    subplot(4,1,2); hold on
    for j = 1:size(Xphi{2},3)
        plot(Xphi{2}(n,:,j,1))
    end
    title('start')
    
    subplot(4,1,3); hold on
    for k = 1:size(Xphi{3},4)
        plot(Xphi{3}(n,:,1,k))
    end
    title('day')
    
    subplot(4,1,4); hold on
    for j = 1:size(Xphi{4},3)
        for k = 1:size(Xphi{4},4)
            plot(Xphi{4}(n,:,j,k))
        end
    end
    title('interaction')
end

