function angle = plot_cov(m, S, varargin)
S = squeeze(S);

[vec, val] = eig(S);
[val, ord] = sort(diag(val));
vec = vec(:,ord);
angle = acos(vec(1,end));
ra = sqrt(val(end));
rb = sqrt(val(1));
ellipse(ra,rb,angle,m(1),m(end), varargin{:});
end