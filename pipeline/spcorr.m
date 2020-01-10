function CORRX = spcorr(X)

[n, k]  = size(X);
Exxprim = full(X'*X)/n; %I'm shocked if this isn't full so let's drop sparse now 
Ex   = full(mean(X))'; %same deal
COVX = (Exxprim - Ex*Ex');
STDEVX = sqrt(diag(COVX));
CORRX = COVX ./ (STDEVX * STDEVX');