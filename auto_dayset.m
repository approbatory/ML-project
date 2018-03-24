function daysets = auto_dayset(names, super)
if nargin == 1
    super = '..';
end
if ~iscell(names)
    names = {names};
end
daysets = cell(1,numel(names));
for ix = 1:numel(names)
    name = names{ix};
    directory = fullfile(super, name);
    S = dir(fullfile(directory,'c*m*d*'));
    for d_ix = 1:numel(S)
        s = S(d_ix);
        [label, short, ds] = ds_autolabel(s.folder, s.name);
        daysets{ix}(d_ix) = struct('directory', s.folder,...
            'day', s.name, 'label', label, 'short', short,...
            'changing', get_changing(short), 'ds', ds, 'res', struct);
    end
end
if numel(names) == 1
    daysets = daysets{1};
end
end

function ch = get_changing(short)
switch short
    case {'RN', 'NR', 'LS', 'SL'}
        ch = 'west';
        return;
    case {'RS', 'SR', 'LN', 'NL'}
        ch = 'east';
        return;
end
ch = '';
end