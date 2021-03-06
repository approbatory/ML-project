function load_from_clouds(o)

basic_vars = {...
    'mus',...
    'evecs',...
    'evecs_shuf',...
    'lambda',...
    'lambda_shuf',...
    'dmus',...
    'loadings',...
    'loadings_shuf',...
    'corr_evecs',...
    'corr_evecs_shuf',...
    'corr_lambda',...
    'corr_lambda_shuf',...
    'corr_loadings',...
    'corr_loadings_shuf'
    };


%progressbar('sess...');
WaitMessage = parfor_wait(o.total_sessions);
for sess_idx = 1:o.total_sessions
    C = Cloud(sess_idx);
    o.mouse{sess_idx} = C.dt.mouse_name;
    
    o.vars.num_neurons{sess_idx} = size(C.dt.data_tensor,1);
    o.vars.num_trials{sess_idx} = size(C.dt.data_tensor,3);
    
    for i = 1:numel(basic_vars)
        varname = basic_vars{i};
        o.vars.(varname){sess_idx} = C.(varname);
        %progressbar([], i/numel(basic_vars));
    end
    
    %progressbar(sess_idx ./ o.total_sessions);
    WaitMessage.Send;
end
WaitMessage.Destroy;

o.save_me;