function [signif_filt, cell_muti, cutoff, prob, shuf_cell_muti] = signif_space_inf(ds, ks, N, given_binary)
K = max(ks);
TAB = tabulate(ks);
counts_y = TAB(:,2);
%Xb = ds_dataset(ds, 'openfield', true, 'filling', 'binary');
%Xb = one_time_events(ds.trials.traces);
Xb = given_binary;
cell_muti = -ones(1,ds.num_cells);
for ci = 1:ds.num_cells
    cell_muti(ci) = muti_bin(Xb(:,ci), ks, K, counts_y);
end

%N = 1000;

shuf_cell_muti = -ones(N,ds.num_cells);
for i = 1:N
    for ci = 1:ds.num_cells
        events_start_shuffled = shuf_x(Xb(:,ci));%es_shuf(ds, ci);
        shuf_cell_muti(i,ci) = muti_bin(events_start_shuffled, ks, K, counts_y);
    end
    %fprintf('i = %d ',i);
    if mod(i,100)==0
        fprintf('i=%d\n',i);
    end
end

cutoff = 0.05 / ds.num_cells;
prob = -ones(1,ds.num_cells);
signif_filt = false(1,ds.num_cells);
for ci = 1:ds.num_cells
    num_over = sum(shuf_cell_muti(:,ci) > cell_muti(ci));
    prob(ci) = num_over/N;
    signif_filt(ci) = prob(ci) < cutoff;
end
end

% function ks_s = shuf_ks(ks)
% ks_s = ks(randperm(length(ks)));
% end

function x_s = shuf_x(x)
evs = sum(x);
tot = length(x);
x_s = sparse(tot,1);
x_s(randperm(tot,evs)) = 1;
end

function es_s = es_shuf(ds, ci)
events = ds.trials.events{ci};
start = events(:,1);
len = events(:,2) - events(:,1);
shuf_start = randperm(ds.full_num_frames, length(start));
%es_s = sparse(ds.full_num_frames, 1);
for ei = 1:length(start)
    if shuf_start(ei) + len(ei) > ds.full_num_frames
        len(ei) = double(ds.full_num_frames) - shuf_start(ei)-1;
    end
end
%one_inds = zeros(sum(len),1); my_place = 1;
ranges_cell = cell(1,length(start));
for ei = 1:length(start)
    %es_s(shuf_start(ei):shuf_start(ei)+len(ei)) = 1;
    ranges_cell{1,ei} = shuf_start(ei):shuf_start(ei)+len(ei);
    %one_inds(my_place:my_place+len(ei)) = shuf_start(ei):shuf_start(ei)+len(ei);
    %one_inds = [one_inds, shuf_start(ei):shuf_start(ei)+len(ei)];
    %my_place = my_place + len(ei) + 1;
end
one_inds = cell2mat(ranges_cell);
if any(one_inds == 0)
    error('bad');
end
es_s = sparse(one_inds, 1, 1, double(ds.full_num_frames), 1);
end