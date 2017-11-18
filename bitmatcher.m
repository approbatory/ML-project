clear all
%profile on
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

[res, outp] = detect_events(traces);

py_tr = load('../ego-fix.tr');

if all(all(res.'==py_tr))
    disp('All match');
else
    disp('not all match');
end

%profile viewer

%%
ix = 123;
rnge = 1:2000;
original = outp{1}(rnge,ix);
submed = outp{2}(rnge,ix);
cleaned = outp{3}(rnge,ix);
valid = outp{4}(rnge,ix);
marked = outp{5}(rnge,ix);

figure
plot(rnge, original);
xlabel frames
title 'original trace'
print original -dpng

figure
plot(rnge, submed);
xlabel frames
title 'median subtracted'
print submed -dpng

figure
plot(rnge, cleaned);
xlabel frames
title 'sliding average'
print cleaned -dpng

figure
plot(rnge, cleaned);
hold on
plot(rnge, valid, 'r');
xlabel frames
title 'valid peaks detected'
print valid -dpng

figure
plot(rnge, cleaned);
hold on
plot(rnge, marked, 'r');
xlabel frames
title 'peaks with width & offset'
print marked -dpng