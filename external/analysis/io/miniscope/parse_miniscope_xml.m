function s = parse_miniscope_xml(source)

s = struct;

% Extract recording parameters
xDoc = xmlread(source);
allAttrs = xDoc.getElementsByTagName('attr');
for k = 0:allAttrs.getLength-1
    attr = allAttrs.item(k);
    attr_type = char(attr.getAttributes.getNamedItem('name').getValue);
    attr_val = char(attr.getFirstChild.getData);
    s = setfield(s, attr_type, attr_val); %#ok<*SFLD>
end

% Extract the list of TIF files. The tag <file> shows up in two 
%   different contexts in the Inscopix XML output. We want only the
%   one that has <decompressed> as its direct parent. Furthermore, strip
%   the absolute path since the files may have moved.
file_names = cell(100,1);
allFiles = xDoc.getElementsByTagName('file');
num_files = 0;
for k = 0:allFiles.getLength-1
    file = allFiles.item(k);
    if strcmp(file.getParentNode.getTagName, 'decompressed')
        num_files = num_files + 1;
        full_path = char(file.getFirstChild.getData);
        [~, name, ext] = fileparts(full_path);
        file_names{num_files} = strcat(name,ext);
    end
end
s = setfield(s, 'num_files', num_files);
s = setfield(s, 'files', file_names(1:num_files));
