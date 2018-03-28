function plot_decoding_results(decoding_results, errtype, errdesc)
if nargin == 1
    errtype = 'dist';
    errdesc = 'bin distance';
end
figure;
for ix = 1:length(decoding_results)
    subplot(1,3,ix);
    label = decoding_results(ix).label;
    res = decoding_results(ix).res;
    algs = cell(1,length(res));
    div_types = cell(1,3);
    for i = 1:length(res)
        algs{i} = res(i).alg.name;
        div = res(i).division;
        for j = 1:length(div)
            div_types{j} = div(j).desc;
            out = div(j).out;
            for k = 1:length(out)
                if strcmp(out(k).err_type, errtype)
                    if iscell(out(k).test_err)
                        d(j,i) = mean(cell2mat(out(k).test_err));
                    else
                        d(j,i) = mean((out(k).test_err));
                    end
                end
            end
        end
    end
    bar(categorical(div_types), mean(d,3));
    title(label);
    ylabel([errdesc ' test']);
    ylim([0 1.5]);
    legend(algs);
end



figure;
for ix = 1:length(decoding_results)
    subplot(1,3,ix);
    label = decoding_results(ix).label;
    res = decoding_results(ix).res;
    algs = cell(1,length(res));
    div_types = cell(1,3);
    for i = 1:length(res)
        algs{i} = res(i).alg.name;
        div = res(i).division;
        for j = 1:length(div)
            div_types{j} = div(j).desc;
            out = div(j).out;
            for k = 1:length(out)
                if strcmp(out(k).err_type, errtype)
                    if iscell(out(k).train_err)
                        d(j,i) = mean(cell2mat(out(k).train_err));
                    else
                        d(j,i) = mean((out(k).train_err));
                    end
                end
            end
        end
    end
    bar(categorical(div_types), mean(d,3));
    title(label);
    ylabel([errdesc ' train']);
    ylim([0 1.5]);
    legend(algs);
end
end