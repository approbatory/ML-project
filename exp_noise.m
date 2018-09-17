ST = SVMTrack('../linear_track/Mouse-2024-20150317_202708-linear-track-TracesAndEvents.mat',...
    14, 20, 50);

%%
%todo make the bins finer

%mean activity in each bin
bin_X = cell(1,ST.num_bins);
for b_ix = 1:ST.num_bins
    bin_X{b_ix} = shuffle(ST.fw_X(ST.fw_ks == b_ix,:), ones(1,sum(ST.fw_ks==b_ix)));
end

mean_bin_X = cellfun(@mean, bin_X, 'UniformOutput', false);
mean_bin_X = cell2mat(mean_bin_X.');

cov_bin_X = cellfun(@cov, bin_X, 'UniformOutput', false);
cov_bin_X = cat(3,cov_bin_X{:});
cov_bin_X = permute(cov_bin_X, [3 1 2]);

y_bin = ST.fw_centers(:,1);
%%

dX = diff(mean_bin_X);
dy = diff(y_bin);
dXdy = dX ./ dy;
princ = zeros(ST.num_bins, ST.total_neurons);
%todo find angle b/w principal axis and dX/dy
angles = zeros(ST.num_bins-1,1);
cos_angles = angles;

for b_ix = 1:ST.num_bins
    Sigma = squeeze(cov_bin_X(b_ix,:,:));
    coeffs = pca(bin_X{1});
    princ(b_ix,:) = coeffs(:,1);
end

for b_ix = 1:ST.num_bins-1
    [angles(b_ix), cos_angles(b_ix)] = angle_v(princ(b_ix,:), dXdy(b_ix,:));
end

angles_tangent = zeros(ST.num_bins-2,1);
cos_angles_tangent = angles;
for b_ix = 1:ST.num_bins-2
    [angles_tangent(b_ix), cos_angles_tangent(b_ix)] = angle_v(dXdy(b_ix,:), dXdy(b_ix+1,:));
end

figure; plot(angles); hold on; plot(angles_tangent)
%%
if false
figure('Position', [100 100 1000 350]);
subplot(1,2,1);
imagesc(squeeze(mpc_c), [0.5 1]);
colorbar;
%axis equal
title 'mean posterior prob. of correct guess'
xlabel bin
ylabel bin
subplot(1,2,2);
imagesc(squeeze(mpc_cs), [0.5 1]);
colorbar;
%axis equal
title 'mean posterior prob. of correct guess - shuffled'
xlabel bin
ylabel bin

m_ = @(C) cellfun(@mean, C);
e_ = @(C) cellfun(@(x) std(x)./sqrt(numel(x)), C);
figure('Position', [1000 1000 1000 350]);
subplot(1,2,1);
errorbar(m_(df_c), e_(df_c));
hold on;
errorbar(m_(df_cs), e_(df_cs));
ylim([0.5 1]);
xlabel 'bin distance'
title 'mean posterior prob. of correct guess (with CV)'
ylabel 'probability'
legend unshuffled shuffled Location best


subplot(1,2,2);
errorbar(m_(df), e_(df));
hold on;
errorbar(m_(df_s), e_(df_s));
ylim([0.5 1]);
xlabel 'bin distance'
title 'mean posterior prob. of correct guess (no CV)'
ylabel 'probability'
legend unshuffled shuffled Location best

end
%%
function [theta, cos_angle] = angle_v(a,b)
a = a(:); b = b(:);
cos_angle = (a.' * b) ./ norm(a) ./ norm(b);
theta = acosd(cos_angle);
end

