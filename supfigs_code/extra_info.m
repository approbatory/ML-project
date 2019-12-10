
cd ..
%%
selected_indices = DecodeTensor.special_by_mouse({'Mouse2022', 'Mouse2024', 'Mouse2028'});
%% calculating regularization
n_reps = 20;
lambda = logspace(-10, -1, 50);

%selected_indices = DecodeTensor.special_by_mouse({'Mouse2022', 'Mouse2024', 'Mouse2028'});

%[accuracy_train, accuracy_test, dp2_train, dp2_test] = ...
%    deal(zeros(n_reps, numel(lambda), numel(selected_indices)));

my_ticker = tic;
parfor i_s_i = 1:numel(selected_indices)
    s_i = selected_indices(i_s_i);
    d = DecodeTensor.cons_filt(s_i);
    
    
    
    for l_i = 1:numel(lambda)
        
        for r_i = 1:n_reps
%             [T1, d1, T2, d2] = DecodeTensor.holdout_half(d.data_tensor, d.tr_dir);
%             
%             [XX1, kk1] = DecodeTensor.tensor2dataset(T1, d1);
%             XX1_sub = XX1(kk1==20 | kk1==22,:);
%             kk1_sub = kk1(kk1==20 | kk1==22,:);
%             
%             [XX2, kk2] = DecodeTensor.tensor2dataset(T2, d2);
%             XX2_sub = XX2(kk2==20 | kk2==22,:);
%             kk2_sub = kk2(kk2==20 | kk2==22,:);
%             
%             alg = my_algs('linsvm', lambda(l_i));
%             
%             model = alg.train(XX1_sub, kk1_sub);
%             pred_train = alg.test(model, XX1_sub);
%             pred_test = alg.test(model, XX2_sub);
%             
%             accuracy_train(r_i, l_i, s_i) = mean(pred_train == kk1_sub);
%             accuracy_test(r_i, l_i, s_i) = mean(pred_test == kk2_sub);
%             
%             dp2_train(r_i, l_i, s_i) = test_model_dp2(model, XX1_sub, kk1_sub);
%             dp2_test(r_i, l_i, s_i) = test_model_dp2(model, XX2_sub, kk2_sub);
            %progressbar(l_i/numel(lambda), r_i/n_reps);
            [accuracy_train{i_s_i}(r_i, l_i), accuracy_test{i_s_i}(r_i, l_i),...
                dp2_train{i_s_i}(r_i, l_i), dp2_test{i_s_i}(r_i, l_i)] = ...
                run_decoding(d, lambda(l_i), []);
        end
        
    end
    
    
    
    
    
    
end
toc(my_ticker)


%% plotting regularization
% figure;
% subplot(2,1,1);
% MultiSessionVisualizer.plot_single_filtered(repmat({lambda},size(accuracy_test)),...
%     {dp2_test, dp2_train}, {'-k', ':k'}, ...
%     [1 2 3]); %selected_indices
% xlim([-Inf Inf]);
% set(gca, 'XScale', 'log');
%%
p = Pub(9, 7, 'rows', 2, 'columns', 2);
p.panel(1, 'xlab', 'L1 regularization \lambda',...
           'ylab', '({\itd}'')^2');
           
line_colors = {'b', 'g', 'k'};
for i_s_i = 1:numel(selected_indices)
    h_1 = Utils.neuseries(lambda, dp2_test{i_s_i}, ['-' line_colors{i_s_i}]);
    %errorbar(lambda, mean(dp2_test{i_s_i}), sem(dp2_test{i_s_i}), ['-' line_colors{i_s_i}]);
    hold on;
    h_2 = Utils.neuseries(lambda, dp2_train{i_s_i}, ['-.' line_colors{i_s_i}]);
    %errorbar(lambda, mean(dp2_train{i_s_i}), sem(dp2_train{i_s_i}), ['-.' line_colors{i_s_i}]);
end
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
%xlabel 'L1 regularization \lambda'
%ylabel '({\itd}'')^2'
ylim([0 Inf]);
legend boxoff
legend([h_1.mainLine h_2.mainLine], 'Test', 'Train', 'Location', 'Best');

p.panel(2, 'xlab', 'L1 regularization \lambda',...
           'ylab', 'Accuracy');
line_colors = {'b', 'g', 'k'};
for i_s_i = 1:numel(selected_indices)
    h_1 = Utils.neuseries(lambda, accuracy_test{i_s_i}, ['-' line_colors{i_s_i}]);
    %errorbar(lambda, mean(accuracy_test{i_s_i}), sem(accuracy_test{i_s_i}), ['-' line_colors{i_s_i}]);
    hold on;
    h_2 = Utils.neuseries(lambda, accuracy_train{i_s_i}, ['-.' line_colors{i_s_i}]);
    %errorbar(lambda, mean(accuracy_train{i_s_i}), sem(accuracy_train{i_s_i}), ['-.' line_colors{i_s_i}]);
end
set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');
%xlabel 'L1 regularization \lambda'
%ylabel Accuracy
ylim([0.5 1]);
legend boxoff
legend([h_1.mainLine h_2.mainLine], 'Test', 'Train', 'Location', 'Best');

p.format;
%%
% subplot(2,1,2);
% errorbar(lambda, mean(accuracy_test(:,:,70)), sem(accuracy_test(:,:,70)), '-', 'DisplayName', 'Test');
% hold on;
% errorbar(lambda, mean(accuracy_train(:,:,70)), sem(accuracy_train(:,:,70)), '-.', 'DisplayName', 'Train');
% set(gca, 'XScale', 'log');
% xlabel 'L1 regularization \lambda'
% ylabel 'Accuracy'
% ylim([0.5 1]);
% legend

%% retrieving dataseries

dbfile = 'decoding_dataseries.db';
conn = sqlite(dbfile);
samp_size = 40;
[sess, mouse_names] = DecodeTensor.filt_sess_id_list;
[data_sizes, imse] = db_imse_reader_datasize(conn, 'unshuffled', sess, samp_size);
[data_sizes_s, imse_s] = db_imse_reader_datasize(conn, 'shuffled', sess, samp_size);
[data_sizes_d, imse_d] = db_imse_reader_datasize(conn, 'diagonal', sess, samp_size);
assert(isequal(data_sizes, data_sizes_s), 'mismatch between unshuffled and shuffled sampling');
assert(isequal(data_sizes, data_sizes_d), 'mismatch between unshuffled and diagonal sampling');

%% plotting dataseries
p.panel(3, 'xlab', 'Number of trials', 'ylab', 'IMSE (cm^{-2})');
MultiSessionVisualizer.plot_single_filtered(data_sizes, {imse_s, imse}, {'r', 'b'}, ...
    selected_indices);
%
%xlabel 'Number of trials'
%ylabel 'IMSE (cm^{-2})'
xlim([-Inf Inf]);
text(80, 0.05, 'Real', 'Color', 'b');
text(80, 0.15, 'Shuffled', 'Color', 'r');
p.format;
%%
for s_i = selected_indices
    d = DecodeTensor.cons_filt(s_i);
    my_ticker = tic;
    parfor rep_i = 1:20
        [sub_data, d1] = DecodeTensor.holdout_half(d.data_tensor, d.tr_dir);
        sub_data_shuf = DecodeTensor.shuffle_tensor(sub_data, d1);
        %sub_data = d.data_tensor;
        [XX_, kk] = DecodeTensor.tensor2dataset(sub_data, d1);
        [XX_s, kk_s] = DecodeTensor.tensor2dataset(sub_data_shuf, d1);
        for n_cl = 1:10
            [~, ~, sumd{rep_i, n_cl}] = kmeans(XX_, n_cl, 'Distance', 'cosine');
            num_samples(rep_i) = size(XX_,1);
            [~, ~, sumd_s{rep_i, n_cl}] = kmeans(XX_s, n_cl, 'Distance', 'cosine');
        end
    end
    toc(my_ticker)
    %
    cos_err{s_i} = 1-cellfun(@sum, sumd)./num_samples';
    cos_err_shuf{s_i} = 1-cellfun(@sum, sumd_s)./num_samples';
end
%cos_err = cellfun(@sum, sumd)./num_samples';
%cos_err_shuf = cellfun(@sum, sumd_s)./num_samples';
%%


%%
p.panel(4, 'xlab', 'Number of clusters, {\itk}',...
    'ylab', sprintf('Mean cosine to\ncluster centroid'));

for s_i = selected_indices
    Utils.neuseries(1:10, cos_err_shuf{s_i}, 'r'); hold on;
    Utils.neuseries(1:10, cos_err{s_i}, 'b'); hold on;
end
%shadedErrorBar(1:10, mean(cos_err_shuf), sem(cos_err_shuf), 'lineprops', 'r'); hold on
%shadedErrorBar(1:10, mean(cos_err), sem(cos_err), 'lineprops', 'b'); hold on;

%xlabel 'Number of clusters, {\itk}'
%ylabel(sprintf('Mean cosine to\ncluster centroid'))
line([2 2], [0 1], 'LineStyle', ':', 'Color', 'k');
ylim([0.25 0.7]);
xlim([1 10]);

%text(6, 0.45, 'Real', 'Color', 'b');
%text(6, 0.35, 'Shuffled', 'Color', 'r');
p.format;
%%
p.print('supplements_pdf', 'extra_info');
%%
function [data_sizes, imse] = db_imse_reader_datasize(conn, setting, sess, samp_size)
%Read from the decoding database.
%Inputs:
%   conn: a valid connection to a sqlite database
%   setting: either 'unshuffled', 'shuffled', or 'diagonal'
%   sess: a string cell of session codes denoting which
%       sessions to read out
%   samp_size: how many samples to load, for a given set of
%       parameters, e.g. 20 or 80
bc = @DecodeTensor.build_command_sess;
q = @Utils.cf_p;
res = q(1,@(s)conn.fetch(bc(s, setting, 'MSE', 'max', [])), sess);
data_sizes = Utils.cf_(@(r)double(cell2mat(r(:,2))), res);
imse = Utils.cf_(@(r)1./cell2mat(r(:,3)), res);
[data_sizes, imse] = cellfun(...
    @(n,i)MultiSessionVisualizer.regroup(n, i, samp_size),...
    data_sizes, imse, 'UniformOutput', false);
end


function dp2 = test_model_dp2(model, XX, kk)
k_values = unique(kk(:));
assert(numel(k_values)==2, 'only binary classification');
decoder_direction = model.Beta;
decoder_direction = decoder_direction(:);
class1 = XX(kk == k_values(1), :);
class2 = XX(kk == k_values(2), :);

class1_proj = class1 * decoder_direction;
class2_proj = class2 * decoder_direction;
assert(isvector(class1_proj) && isvector(class2_proj), 'not vectors!');
dp2 = (mean(class1_proj) - mean(class2_proj)).^2 ./...
    ((var(class1_proj) + var(class2_proj))./2);
end

function [accuracy_train, accuracy_test, dp2_train, dp2_test] = ...
    run_decoding(d, lam, data_size)
[T, d] = DecodeTensor.cut_tensor(d.data_tensor, d.tr_dir, [], data_size);
[T1, d1, T2, d2] = DecodeTensor.holdout_half(T, d);

[XX1, kk1] = DecodeTensor.tensor2dataset(T1, d1);
XX1_sub = XX1(kk1==20 | kk1==22,:);
kk1_sub = kk1(kk1==20 | kk1==22,:);

[XX2, kk2] = DecodeTensor.tensor2dataset(T2, d2);
XX2_sub = XX2(kk2==20 | kk2==22,:);
kk2_sub = kk2(kk2==20 | kk2==22,:);

alg = my_algs('linsvm', lam);

model = alg.train(XX1_sub, kk1_sub);
pred_train = alg.test(model, XX1_sub);
pred_test = alg.test(model, XX2_sub);

accuracy_train = mean(pred_train == kk1_sub);
accuracy_test = mean(pred_test == kk2_sub);

dp2_train = test_model_dp2(model, XX1_sub, kk1_sub);
dp2_test = test_model_dp2(model, XX2_sub, kk2_sub);
end