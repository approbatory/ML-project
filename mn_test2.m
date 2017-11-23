function pred = mn_test2(model,X_test)
[~,pred] = max(model.log_prior + X_test*model.log_conditional);
end

