function load_definitions(org)

org.make_derived('signal_variance', {'lambda', 'loadings'},...
    @(lambda, loadings) lambda .* loadings'.^2);

org.make_derived('corr_signal_variance', {'corr_lambda', 'corr_loadings'},...
    @(lambda, loadings) lambda .* loadings'.^2);

org.make_derived('signal_density', {'dmus'},...
    @(dmu) (sum(dmu.^2).^2 ./ sum(dmu.^4)) ./ length(dmu));

org.make_derived('cos_area', {'loadings', 'loadings_shuf'},...
    @area_between);

org.make_derived('cos2_area', {'loadings', 'loadings_shuf'},...
    @(l,ls)area_between(l.^2, ls.^2));

org.make_derived('cos_area_per_neuron', {'cos_area', 'num_neurons'},...
    @(a,b)a./b);

org.make_derived('cos2_area_per_neuron', {'cos2_area', 'num_neurons'},...
    @(a,b)a./b);


org.make_derived('intersect', {'loadings', 'loadings_shuf'},...
    @(a,b)find(a<b,1));
end


function r = area_between(a,b)

%int = find(a < b, 1);
int = [];

if isempty(int)
    r = sum(abs(a-b));
else
    r = sum(abs(a(1:int-1) - b(1:int-1)));
end
end