function s = fastdlda

s.name = 'Fast DLDA';
s.pre = @(x)x;
s.train = @my_dlda;
s.test = @my_test_dlda;