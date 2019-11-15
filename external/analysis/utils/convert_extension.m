function new_filestr = convert_extension(filestr, new_ext)

[pathstr, name, ~] = fileparts(filestr);
new_filestr = fullfile(pathstr, strcat(name, '.', new_ext));