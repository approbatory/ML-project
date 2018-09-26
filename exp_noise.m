ST = SVMTrack('../linear_track/Mouse-2024-20150317_202708-linear-track-TracesAndEvents.mat',...
    14, 20, 20);
%%
% basic
[bin_X, y_bin, mean_bin_X, cov_bin_X] = SVMTrack.bin_data(ST.fw_X, ST.fw_ks, ST.fw_centers, false, true);
[bin_X_shuf, y_bin_shuf, mean_bin_X_shuf, cov_bin_X_shuf] = SVMTrack.bin_data(ST.fw_X, ST.fw_ks, ST.fw_centers, true, true);

% with PLS
[XL, ~, X_pls] = plsregress(ST.fw_X, ST.fw_scy(:,1), 2);
[bin_X_pls, y_bin_pls, mean_bin_X_pls, cov_bin_X_pls] = SVMTrack.bin_data(X_pls, ST.fw_ks, ST.fw_centers, false, true);

X_shuf = shuffle(ST.fw_X, ST.fw_ks);
[~, ~, X_shuf_pls] = plsregress(X_shuf, ST.fw_scy(:,1), 2);
[bin_X_shuf_pls, y_bin_shuf_pls, mean_bin_X_shuf_pls, cov_bin_X_shuf_pls] = SVMTrack.bin_data(X_shuf_pls, ST.fw_ks, ST.fw_centers, false, true);

%%
colors = parula(1000);
c_ix = 1+round(999.*(y_bin_pls - min(y_bin_pls))./range(y_bin_pls));

figure;
subplot(2,2,1);
hold on
scatter(X_pls(:,1), X_pls(:,2), 1, ST.fw_scy(:,1));
scatter(mean_bin_X_pls(:,1), mean_bin_X_pls(:,2), 10, [1 0 0]);%y_bin_pls);
xlim_ = xlim;
ylim_ = ylim;
xlim(xlim_);
ylim(ylim_);
xlabel PLS1
ylabel PLS2
title(['Linear track activity, ' num2str(ST.num_bins) ' bins']);

subplot(2,2,2);
hold on;
for b_ix = 1:ST.num_bins
    plot_cov(mean_bin_X_pls(b_ix,:), cov_bin_X_pls(b_ix,:,:), colors(c_ix(b_ix),:));
end
xlim(xlim_);
ylim(ylim_);
xlabel PLS1
ylabel PLS2
title(['Visualized covariances, ' num2str(ST.num_bins) ' bins']);

subplot(2,2,3);
hold on
scatter(X_shuf_pls(:,1), X_shuf_pls(:,2), 1, ST.fw_scy(:,1));
scatter(mean_bin_X_shuf_pls(:,1), mean_bin_X_shuf_pls(:,2), 10, [1 0 0]);%y_bin_pls);
xlim(xlim_);
ylim(ylim_);
xlabel PLS1
ylabel PLS2
title(['Linear track activity - shuffled, ' num2str(ST.num_bins) ' bins']);

subplot(2,2,4);
hold on;
for b_ix = 1:ST.num_bins
    plot_cov(mean_bin_X_shuf_pls(b_ix,:), cov_bin_X_shuf_pls(b_ix,:,:), colors(c_ix(b_ix),:));
end
xlim(xlim_);
ylim(ylim_);
xlabel PLS1
ylabel PLS2
title(['Visualized covariances - shuffled, ' num2str(ST.num_bins) ' bins']);
%%

% dX = diff(mean_bin_X);
% dy = diff(y_bin);
% dXdy = dX ./ dy;
% princ = zeros(ST.num_bins, ST.total_neurons);
% %todo find angle b/w principal axis and dX/dy
% angles = zeros(ST.num_bins-1,1);
% cos_angles = angles;
% 
% for b_ix = 1:ST.num_bins
%     Sigma = squeeze(cov_bin_X(b_ix,:,:));
%     coeffs = pca(bin_X{1});
%     princ(b_ix,:) = coeffs(:,1);
% end
% 
% for b_ix = 1:ST.num_bins-1
%     [angles(b_ix), cos_angles(b_ix)] = angle_v(princ(b_ix,:), dXdy(b_ix,:));
% end
% 
% angles_tangent = zeros(ST.num_bins-2,1);
% cos_angles_tangent = angles;
% for b_ix = 1:ST.num_bins-2
%     [angles_tangent(b_ix), cos_angles_tangent(b_ix)] = angle_v(dXdy(b_ix,:), dXdy(b_ix+1,:));
% end
[angles, angles_tangent] = ST.find_angles(mean_bin_X, bin_X, y_bin);
[angles_shuf, angles_tangent_shuf] = ST.find_angles(mean_bin_X_shuf, bin_X_shuf, y_bin_shuf);
%%

figure;
hold on;

plot(angles, '-o');
plot(angles_shuf, '-o');
plot(angles_tangent, '-o');
L_ = refline(0, 90);
L_.Color = 'k';
%plot(angles_tangent_shuf, 'r.-')

legend 'principal noise vs. tangent, unshuffled' 'principal noise vs. tangent, shuffled' '\Delta tangent angle' 'orthogonal'
xlabel 'bin'
ylabel 'angle (degrees)'
title 'Angle between maximal noise direction and direction of change in position parameter';
%% using bin_X and mean_bin_X, do the rotation scheme

dX = diff(mean_bin_X);

res_bin_X = cell(1,20);
rot_res_bin_X = cell(1,20);
rotator_bin = cell(1,20);
for i = 1:20
    res_bin_X{i} = bin_X_shuf{i} - mean_bin_X(i,:);
    rotator_bin{i} = rotator(dX(1,:).', dX(min(i,19),:).');
    rot_res_bin_X{i} = res_bin_X{i} * rotator_bin{i};
end

pooled_residuals = cell2mat(rot_res_bin_X.');
coeffs = pca(pooled_residuals);
princ = coeffs(:,1);
[theta, cos_theta] = angle_v(princ, dX(1,:));
fprintf(['The angle between principal noise component (pooled over 20 bins)'...
    ', and the direction of increase is:\n%.2f degrees\n'], theta);
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
plot_cov_ = @(m,S) plot_cov(m,S);
%%
function [theta, cos_angle] = angle_v(a,b)
a = a(:); b = b(:);
cos_angle = (a.' * b) ./ norm(a) ./ norm(b);
theta = acosd(cos_angle);
end

function plot_cov(m, S, varargin)
S = squeeze(S);

[vec, val] = eig(S);
[val, ord] = sort(diag(val));
vec = vec(:,ord);
angle = acos(vec(1,end));
ra = sqrt(val(end));
rb = sqrt(val(1));
ellipse(ra,rb,angle,m(1),m(end), varargin{:});
end