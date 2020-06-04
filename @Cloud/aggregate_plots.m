function [N_50_dist, sig_agg, noise_agg, signoise_agg, source_mouse] = aggregate_plots
N = 107;
nb = 40;
N_50_dist = zeros(1,N);
sig_agg = zeros(nb-1, nb-1, N);
noise_agg = zeros(nb, nb, N);
signoise_agg = zeros(nb, nb-1, N);

source_mouse = cell(1,N);

valid = false(1,N);
progressbar('agg plots');
for i = 1:N
    %try
        C = Cloud(i);
        N_50_dist(i) = C.N;
        sig_agg(:,:,i) = C.signal_geo;
        noise_agg(:,:,i) = C.noise_geo;
        signoise_agg(:,:,i) = C.signal_noise_overlap_geo;
        source_mouse{i} = C.dt.mouse_name;
        valid(i) = true;
        fprintf('Index %d succeeded\n', i);
    %catch
    %    fprintf('index %d failed\n', i);
    %end
    progressbar(i/N);
end
N_50_dist = N_50_dist(valid);
sig_agg = sig_agg(:,:,valid);
noise_agg = noise_agg(:,:,valid);
signoise_agg = signoise_agg(:,:,valid);
source_mouse = source_mouse(valid);

figure;
%scatter(categorical(source_mouse), N_50_dist);
boxplot(N_50_dist, source_mouse);
xlabel 'Mouse'
ylabel 'Chosen {\itn}_{1/2} in each session'

figure;
subplot(1,4,1);
histogram(N_50_dist);
xlabel 'Chosen subspace dim'
title(sprintf('Noise subspace dimension\ncontaining 50%% of \\Delta\\mu'));

subplot(1,4,2);
imagesc(mean(sig_agg,3));
xlabel 'Spatial bin'
ylabel 'Spatial bin'
colormap(gca, bluewhitered);
colorbar;
title 'Cos overlap between signal directions'
axis image

subplot(1,4,3);
imagesc(mean(noise_agg,3));
line([20 20]+0.5, ylim, 'Color', 'w');
line(xlim, [20 20]+0.5, 'Color', 'w');
xlabel 'Spatial bin';
ylabel 'Spatial bin';
colorbar;
title(sprintf('Mean canon corr between\nnoise subspaces'));
axis image

subplot(1,4,4);
imagesc(mean(signoise_agg,3));
line([20 20]+0.5, ylim, 'Color', 'w');
line(xlim, [20 20]+0.5, 'Color', 'w');
xlabel 'Spatial bin of signal direction';
ylabel 'Spatial bin of noise subspace';
colorbar;
title(sprintf('Cos overlap between noise subspaces\nand signal directions'))
axis image;
end