% h = FunctionHandleFromPath(pathToFunction)
%
% Takes the path to an m-file containing a function definition and
% returns a handle to the function.
function hf = FunctionHandleFromPath(pathToFunction)
clear functions;
currentDir = pwd;
[pth, name] = fileparts(pathToFunction);
cd(pth);
hf = str2func(name);
cd(currentDir);