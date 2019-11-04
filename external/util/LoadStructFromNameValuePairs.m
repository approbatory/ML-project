function S = LoadStructFromNameValuePairs(cellArray, varargin)
% S = LoadStructFromNameValuePairs(cellArray, varargin)
%
% Loads the parameters specified as name/value pairs in 'cellArray'
% into a structure. An optional list of fieldnames can be provided,
% which will be filled from the data in the cell array. If the
% required fields are missing, an error will be generated, unless
% a second cell array of default values is provided.
%
% Example:
%
% S = LoadStructFromNameValuePairs({'firstName','Michael','lastName','Smith','age', 31})
%
% After this call, S looks like:
% S.firstName = 'Michael'
% S.lastName = 'Smith';
% S.age = 31;
%

S = struct;
if (mod(numel(cellArray),2)~=0)
  error('Input argument must have an even number of elements.');
end

for i = 1:2:numel(cellArray)
  fieldName = cellArray{i};
  fieldValue= cellArray{i+1};
  if (~ischar(fieldName))
    error('Element %d of the cell array is not a string.', i);
  end
  S.(fieldName) = fieldValue;
end

switch(numel(varargin))
  case 1
   for i = 1:numel(varargin)
     if (~isfield(S, varargin{i}))
       error('Could not find data for required field "%s".', varargin{i});
     end
   end
 case 2
  fields = varargin{1};
  values = varargin{2};
  for i = 1:numel(fields)
    if (~isfield(S, fields{i}))
      S.(fields{i}) = values{i};
    end
  end
end


    
