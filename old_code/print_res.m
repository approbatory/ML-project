function print_res(res, errtype)
if nargin == 1
    errtype = 'dist';
end
for i = 1:length(res)
    fprintf('%s\n',res(i).alg.name);
    for j = 1:length(res(i).division)
        fprintf('\t%s\t',res(i).division(j).desc);
        for k = 1:length(res(i).division(j).out)
            if strcmp(res(i).division(j).out(k).err_type, errtype)
                fprintf('train_err: %f\ttest_err: %f',...
                    res(i).division(j).out(k).train_err,...
                    res(i).division(j).out(k).test_err);
            end
        end
        fprintf('\n');
    end
end
end