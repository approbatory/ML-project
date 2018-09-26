function ret = rotator(a, b)
assert(isvector(a) && isvector(b));
assert((size(a,2) == 1) && (size(b,2) == 1)); 
A = [a,b];
[Q, R] = qr(A);

scaa = R(1,1);
vecb = R(1:2,2);

rotnum = vecb./scaa;
rotnum = rotnum ./ norm(rotnum);

c = rotnum(1); s = rotnum(2);

rot = [c s; -s c];

Nt = [Q(:,1:2) * rot, Q(:,3:end)];
ret = Q*Nt.';
end