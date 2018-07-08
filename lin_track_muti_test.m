function [signif_filt, cell_muti] = lin_track_muti_test(ds, N, openf)
if ~exist('N', 'var')
    N = 10000;
end

if ~exist('openf', 'var')
    openf = false;
end


%X = ds.trials.traces(:,sel_fw).';
if openf
    y = ds.trials.centroids;
    ks = gen_place_bins(y, 8, 46);
    X_rd = tevents(ds.trials.traces);
    given_binary = X_rd;
    [signif_filt, cell_muti] = signif_space_inf(ds, ks, N, given_binary);

else
    [sel_fw, sel_bw] = select_directions(ds.trials.centroids);
    y = ds.trials.centroids(sel_fw,:);
    ks = gen_place_bins(y, 20, 120*0.8);
    X_rd = tevents(ds.trials.traces);
    given_binary = X_rd(sel_fw,:);
    [signif_filt_fw, cell_muti_fw] = signif_space_inf(ds, ks, N, given_binary);

    y = ds.trials.centroids(sel_bw,:);
    %X = ds.trials.traces(:,sel_bw).';
    ks = gen_place_bins(y, 20, 120*0.8);
    X_rd = tevents(ds.trials.traces);
    given_binary = X_rd(sel_bw,:);
    [signif_filt_bw, cell_muti_bw] = signif_space_inf(ds, ks, N, given_binary);

    signif_filt = [signif_filt_fw ; signif_filt_bw];
    cell_muti = [cell_muti_fw ; cell_muti_bw];
end

end

