function [stabXY, XY] = stabilize_rotor(ds, pro)
XY = cell2mat(preprocess_xy(ds));
X = XY(:,1); Y = XY(:,2);
%horiz = XY(abs(X) > 2*abs(Y),2); %try to get the y values of horiz. points to be close to 0
%vert = XY(abs(Y) > 2*abs(X),1); %try to get the x values of vert. points to be close to 0

%data = [horiz, ones(size(horiz)); vert, ones(size(vert))];
% [~, ord] = sort(abs(X) - abs(Y));
% sel = floor(length(ord)*pro);
% top = XY(ord(1:sel),:);
% bottom = XY(ord(end-sel:end),:);
% data = [top; bottom];


East = XY(X > (1-pro),:);
West = XY(X < pro,:);
North = XY(Y > (1-pro),:);
South = XY(Y < pro,:);

data = [East; West; North; South] - 0.5;
%coeffs = pca(data);
%if det(coeffs) < 0
%    coeffs = [0 1;1 0] * coeffs;
%end
%stabXY = (XY - 0.5) * coeffs + 0.5;
res = rica(data,2);
if det(res.TransformWeights) < 0
    res.TransformWeights = [0 1;1 0] * res.TransformWeights;
end
stabXY = res.transform(XY - 0.5) + 0.5;
end