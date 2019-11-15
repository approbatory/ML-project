function UnpackStructureFieldsAsVariables(S)
% UnpackStructureFieldsAsVariables(S)
%
% For each field in S, creartes a variable in the caller workspace
% with the name of the field and sets its value to that of the field.

fieldNames = fieldnames(S);

for i = 1:numel(fieldNames)
  assignin('caller', fieldNames{i}, S.(fieldNames{i}));
end
