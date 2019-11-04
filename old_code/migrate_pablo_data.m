%starting from the pablo's data, migrate it to the ds
%framework, so that data can be processed uniformly

S = dir('../open_field/_mouse2028_full');
ix = 0;
for i = 1:numel(S)
    if ~S(i).isdir
        ix = ix+1;
        ticker = tic;
        E_T{ix} = load(fullfile(S(i).folder, S(i).name));
        toc(ticker);
        fprintf('Loaded day %d\n', ix);
    end
end


full_num_frames = cell(1,numel(E_T));
num_cells = cell(1,numel(E_T));
trials = cell(1,numel(E_T));
for ix = 1:numel(E_T)
    y = E_T{ix}.tracesEvents.position; %(T, 2)
    X = E_T{ix}.tracesEvents.rawTraces; %(T, C)
    if ix == 1
        y = y(1:6000,:);
        X = X(1:6000,:);
    elseif ix == 4
        y = y(1000:end,:);
        X = X(1000:end,:);
    elseif ix == 12
        y = y(3376:end,:);
        X = X(3376:end,:);
    end
    
    full_num_frames{ix} = size(y,1);
    num_cells{ix} = size(X,2);
    trials{ix}.centroids = y;
    trials{ix}.traces = X.';    
end

my_ds = struct('full_num_frames', full_num_frames,...
    'num_cells', num_cells, 'trials', trials);
fprintf('Created struct array, saving...\n');
save '../linear_track/pablo_data_ds.mat' my_ds
fprintf('Done\n');