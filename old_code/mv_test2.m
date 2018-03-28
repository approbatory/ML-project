function pred = mv_test2(model, X_test)
%M = size(X_test,1);
%pred = zeros(M,1);
%for i = 1:M
%    X = X_test(i,:);
%    l = model.log_prior' + sum(model.log_conditional(:, X==1),2) + sum(log1p(-exp(model.log_conditional(:,X==0))),2);
%    [~,pred(i)] = max(l);
%end
[~, pred] = max(model.log_prior' + sum(model.log_flip_conditional,2) +...
    (model.log_conditional - model.log_flip_conditional)*X_test');
pred = pred';
pred = model.K_vals(pred);
end