function args = UnpackStructureFieldsAsNameValuePairs(S)
% args = UnpackStructureFieldsAsNameValuePairs(S)
%
% Unpacks the fields and values of S into consecutive elements of a
% cell array so that they can be passed as Name Value pairs to
% functions. 
%
% Example:
%
% >> S = struct;
% >> S.fname = 'John';
% >> S.lname = 'Smith';
% >> UnpackStructureFieldsAsNameValuePairs(S)
%
% ans = 
%
%     'fname'    'John'    'lname'    'Smith'
%
% >> whos ans
%   Name      Size            Bytes  Class    Attributes
%
%   ans       1x4               486  cell               


fields = fieldnames(S);
args = {};
for i = 1:numel(fields)
  args{2*i-1} = fields{i};
  args{2*i}   = S.(fields{i});
end

