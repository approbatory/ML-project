function var = LoadVarFromMatFileByName(matFile, varName)
% var = LoadVarFromMatFileByName(matFile, varName)
%
% Loads the variable named 'varName' from 'matFile'. 

data = load(matFile,varName);
var = data.(varName);
