function view_errmap(ds, poss, err_maps, vals, masks, label, alg_label)
class_to_label = containers.Map({1,2,3,4},{'north','south','east','west'});
POS_LABEL = 'position';
num_partitions = length(err_maps);
correct_trials = {ds.trials.correct};
trials = cell(1,num_partitions);
for k = 1:num_partitions
    trials{k} = correct_trials(masks{k});
end
figure;
for k = 1:num_partitions
    subplot(1,num_partitions,k);
    err_map = err_maps{k};
    imagesc(1 - err_map - (1-err_map).*0.1.*mod((1:size(err_map,1))',2), [0 1]);
    colormap('gray');
    [n_trials, n_pos] = size(err_map);
    for i = 1:n_trials
        if trials{k}{i}
            c = 'g';
        else
            c = 'r';
        end
        rectangle('Position', [(n_pos-0.5) (i-0.5) 2 1], 'FaceColor', c);
    end
    ticks = get(gca, 'XTick');
    set(gca, 'XTick', ticks, 'XTickLabel', poss(ticks));
    %t = text(n_pos+2, n_trials/4, 'non-rewarded in red');
    %t.Rotation = -90;
    xlabel(POS_LABEL);
    ylabel('Trials');
    title(sprintf('%s | %s: %s',class_to_label(vals(k)), label, alg_label));
end

end