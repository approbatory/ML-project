t_ = tic;
max_rate = 0.25;
g_field = @(x, m, s) max_rate.*exp(-((x-m)./s).^2 ./ 2);

n_cells = 50e3;
track_len = 120;
field_std = 3;

x = 0:0.1:track_len; % cm
m = sort(rand(n_cells, 1) * track_len);

%R = g_field(x, m, field_std);

rate_func = @(x) g_field(single(x), single(m), single(field_std));
response_sample = @(x) poissrnd(rate_func(x));

t = seconds(seconds(0):seconds(1/20):hours(0.5)); %s
max_v = 30; %cm/s
x_traj = track_len .* sin(max_v./track_len .* t).^2;

noise_sigma = 10;
input_noise = noise_sigma*randn(size(x_traj));

X_noise = uint8(response_sample(x_traj + input_noise));

tracesEvents.simulated_noise = X_noise.';
tracesEvents.position = [x_traj(:), 0*x_traj(:)];
savefile = 'Mouse-2099-20210606_000000-linear-track-TracesAndEvents.mat';
save(savefile, 'tracesEvents');


d_sim_noise = DecodeTensor({savefile, 'Mouse2099'}, 'simulated_noise');

%%
n_cells = size(d_sim_noise.data_tensor, 1);

n_neu = unique(round(logspace(0,log10(n_cells),50)));

progressbar('reps...', 'neu...');
n_reps = 20;

[mse, mse_sh] = deal(zeros(n_reps, numel(n_neu)));

for r_ix = 1:n_reps
for n_ix = 1:numel(n_neu)
    
    err_res = DecodeTensor.decode_all(d_sim_noise.data_tensor,...
        d_sim_noise.tr_dir, 6, fastpnb, n_neu(n_ix), []); 
    
    mse(r_ix, n_ix) = err_res.MSE.unshuffled;
    mse_sh(r_ix, n_ix) = err_res.MSE.shuffled;
    
    progressbar([],n_ix/numel(n_neu));
end

progressbar(r_ix/n_reps);
end

%
fprintf('Doing short...\n');
n_cells = 500;

d_neu = 10;
n_neu_short = unique([1, d_neu:d_neu:n_cells, n_cells]);

progressbar('reps...', 'neu...');
n_reps = 20;

[mse_short, mse_sh_short] = deal(zeros(n_reps, numel(n_neu_short)));


for r_ix = 1:n_reps
for n_ix = 1:numel(n_neu_short)
    err_res = DecodeTensor.decode_all(d_sim_noise.data_tensor,...
        d_sim_noise.tr_dir, 6, fastpnb, n_neu_short(n_ix), []);
    
    mse_short(r_ix, n_ix) = err_res.MSE.unshuffled;
    mse_sh_short(r_ix, n_ix) = err_res.MSE.shuffled;
    
    progressbar([],n_ix/numel(n_neu_short));
end

progressbar(r_ix/n_reps);
end

%%
fr_real = createFit_infoSaturation(n_neu, mean(1./mse));
fr_real_short = createFit_infoSaturation(n_neu_short, mean(1./mse_short));

fr_shuf = createFit_infoSaturation(n_neu, mean(1./mse_sh));
fr_shuf_short = createFit_infoSaturation(n_neu_short, mean(1./mse_sh_short));

figure('Position', [472 626 1048 471]);
subplot(2,4,[1 5]);
hr = ebI(1, fr_real, 'b'); hold on;
ebI(2, fr_real_short, 'b');
hs = ebI(3, fr_shuf, 'r'); hold on;
ebI(4, fr_shuf_short, 'r');
%set(gca, 'YScale', 'log');
xlim([0 5]);
%ylim([1 Inf]);
ylabel '{\itI}_0 fit value (cm^{-2}/neuron)'
set(gca, 'XTick', [1 2 3 4]);
set(gca, 'XTickLabels', {'50k', '500', ...
    '50k shuf', '500 shuf'});
% legend([hr hs], {'Real' 'Shuffled'});
% legend boxoff
% legend location northwest
Utils.fix_exponent(gca, 'y', 1);

subplot(2,4,[2 6]);
hr = eb(1, fr_real, 'b'); hold on;
eb(2, fr_real_short, 'b');
hs = eb(3, fr_shuf, 'r'); hold on;
eb(4, fr_shuf_short, 'r');
set(gca, 'YScale', 'log');
xlim([0 5]);
ylim([1 Inf]);
ylabel '{\itN} fit value (neurons)'
set(gca, 'XTick', [1 2 3 4]);
set(gca, 'XTickLabels', {'50k', '500', ...
    '50k shuf', '500 shuf'});
% legend([hr hs], {'Real' 'Shuffled'});
% legend boxoff
% legend location northwest

%
subplot(2,4,[3 4]);
ep(n_neu, 1./mse_sh, 'r.'); hold on;
plot(n_neu, fr_shuf(n_neu), 'r');

ep(n_neu, 1./mse, 'b.'); 
plot(n_neu, fr_real(n_neu), 'b');

xlabel 'Num neurons'
ylabel '1/MSE (cm^{-2})'
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title 'Up to 50,000 neurons'

text(2, 20, 'Simulated', 'Color', 'b');
text(2, 5, 'Simulated & shuffled', 'Color', 'r');


subplot(2,4,[7 8]);
ep(n_neu_short, 1./mse_sh_short, 'r.'); hold on;
plot(n_neu_short, fr_shuf_short(n_neu_short), 'r');

ep(n_neu_short, 1./mse_short, 'b.');
plot(n_neu_short, fr_real_short(n_neu_short), 'b');
xlabel 'Num neurons'
ylabel '1/MSE (cm^{-2})'
title 'Up to 500 neurons'

sgtitle 'Parameter fits on 50,000 vs. 500 simulated neurons'

function b = eb(ix, sf, c)
    ci = confint(sf);
    ci = ci(:,2);
    e = num2cell(ci - sf.N);
    b = bar(ix, sf.N, c); hold on;
    errorbar(ix, sf.N, e{:}, 'k');
end

function b = ebI(ix, sf, c)
    ci = confint(sf);
    ci = ci(:,1);
    e = num2cell(ci - sf.I_0);
    b = bar(ix, sf.I_0, c); hold on;
    errorbar(ix, sf.I_0, e{:}, 'k');
end

function ep(n, m, varargin)
m_mean = mean(m, 1);
samples = size(m, 1);
m_std = std(m, [], 1);
m_sem = m_std ./ sqrt(samples);
errorbar(n, m_mean, m_sem.*1.96, varargin{:},...
    "CapSize", 1, "Marker", "none", 'LineStyle', "none");
end