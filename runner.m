load ~/assets/c14m4d15/cm01-fix/rec_161217-182928.mat
lines = strsplit(fileread('/home/omer/assets/c14m4d15/cm01-fix/class_161217-183059.txt'), '\n');
goodcells = [];
for l=lines
    if ~isempty(l{1}) && isempty(strfind(l{1}, 'not'))
        goodcells = [goodcells, sscanf(l{1}, '%d')];
    end
end

traces = traces(:, goodcells);


%fprintf('sum of 123 is %f\n', sum(traces(:,ix)));

[res, outp] = detect_events(traces.');
%%
%5852 and 8080 are omitted here
display(find(outp{4}(310,:)))
figure
plot(outp{3}(310,:), '-x')
hold on
plot(std(outp{3}(310,:)) * ones(1,length(outp{3}(310,:))))
%%
tr = outp{3}(310,:);
ix = 12142;
display(tr(ix));
thresh = 1.5*std(tr);
display(thresh);
passing_thresh = tr >= thresh;
display(passing_thresh(ix));
repeated_passing = conv2(1.0*passing_thresh, 1.0*[1,1,1,1,1], 'same') >= 5;
display(repeated_passing(ix));
local_maximum = [1, diff(diff(tr, 1, 2) > 0, 1, 2) < 0, 1];
display(local_maximum(ix));
well_sep_max = local_maximum;
for i = 1:4
    %well_sep_max = well_sep_max & ~circshift(local_maximum, [0,i]);
    well_sep_max = well_sep_max & ~circshift(well_sep_max, [0,i]);
end
display(well_sep_max(ix));
[retest, nuo] = detect_events(outp{1}(310,:));
%display(retest(ix-10:ix));
display(nuo{4}(ix));
%%
for ix = 1:length(goodcells)
fprintf('%d transients sum: %f\n', ix, (sum(outp{1}(ix,:))));
fprintf('%d sub med sum: %f\n', ix, (sum(outp{2}(ix,:))));
fprintf('%d transients_cleaned sum: %f\n', ix, (sum(outp{3}(ix,:))));
fprintf('%d transients_maxima sum: %d\n', ix, int64(sum(outp{4}(ix,:))));
%found_peaks = find(outp{4}(ix,:));
%display(found_peaks);
fprintf('%d transients_canvas sum: %d\n', ix, sum(outp{5}(ix,:)));
fprintf('\n\n');
end
%%
fprintf('%d transients sum: %f\n', ix, sum(sum(outp{1})));
fprintf('%d sub med sum: %f\n', ix, sum(sum(outp{2})));
fprintf('%d transients_cleaned sum: %f\n', ix, sum(sum(outp{3})));
fprintf('%d transients_maxima sum: %d\n', ix, int64(sum(sum(outp{4}))));
%found_peaks = find(outp{4}(ix,:));
%display(found_peaks);
fprintf('%d transients_canvas sum: %d\n', ix, sum(sum(outp{5})));
fprintf('\n\n');
%%
[rowz, colz] = find(outp{4});

[sorted_rowz, indices] = sort(rowz);

combined = [rowz(indices), colz(indices)];

%checking m not present in p
[num_rows_m, num_cols_m] = size(combined);
[num_rows_p, num_cols_p] = size(py_max);
 
display(setxor(combined, py_max, 'rows'));

display(ismember([310 12142], py_max, 'rows'));
display(ismember([310 12142], combined, 'rows'));

%%
display('items in m but not in p');
for i = 1:num_rows_m
    entry = combined(i,:);
    present = 0;
    for j = 1:num_rows_p
        if all(entry==py_max(j,:))
            present = 1;
        end
    end
    if ~present
        display(entry);
    end
end