function s = esc(s)
s = replace(s, '\', '\\');
s = replace(s, '_', '\_');
s = replace(s, '^', '\^');
end