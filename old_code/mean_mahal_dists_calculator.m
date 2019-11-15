function mean_mahal_dists = mean_mahal_dists_calculator(models_ord)
%calculate average mahal
num_folds = length(models_ord{1});
num_bins = size(models_ord{1}{1}.Mu,1);
mean_mahal_dists = zeros(1,num_folds);
for f_ix = 1:num_folds
    my_model = models_ord{1}{f_ix};
    sum_mahal_dists = 0;
    count_mahal_dists = 0;
    for b_ix = 1:num_bins
        mahal_dists = my_model.mahal(my_model.Mu(b_ix,:));
        sum_mahal_dists = sum_mahal_dists + sum(mahal_dists);
        count_mahal_dists = count_mahal_dists + (num_bins-1);
    end
    mean_mahal_dists(f_ix) = sum_mahal_dists/count_mahal_dists;
end