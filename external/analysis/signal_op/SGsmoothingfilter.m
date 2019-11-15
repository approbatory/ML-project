function h = SGsmoothingfilter(N,L)
% function [h,SGfilter] = SGsmoothingfilter(N,L);
%
% Design Savitzky - Golay smoothing filter 
%  with (2N+1) points and L'th order
%
% (Projects (2N+1) dimensional input to 
%  the L'th order polynomial space)
% 
% Output : 
%    h : Impulse response of the filter
%    SGfilter : Expression of the filter in terms of finite differences
%
% August 2013, 
% Cagatay Candan
%

if L>(2*N), disp('Error : L <= 2*N should be satisfied'); return; end;

p = zeros(2*N+1,1); p(N+1)=1;  % particular solution for smoothing app.
%syms D; 
dimnull = 2*N-L;

vec=1; Nullmat = zeros(2*N+1,2*N);
Nullmat(N,1)=1/2;Nullmat(N+2,1)=-1/2;
for dum=2:2:2*N
   vec = conv([1 -2 1],vec); 
   coln = Nullmat(:,dum);
   coln(N+1-dum/2:N+1+dum/2) = vec;
   Nullmat(:,dum) = coln;
   if dum<2*N,
       coln = Nullmat(:,dum+1);
       coln(N-dum/2:N+2+dum/2) = conv(vec,[-1/2 0 1/2]);
       Nullmat(:,dum+1) = coln;
   end;
end;

Nullmat = Nullmat(:,(L+1):end);
Nullmat = sym(Nullmat);

if ~isempty(Nullmat)    
    coefs = (Nullmat'*Nullmat)\Nullmat'*p;
    %SGfilter = 1 - (D*ones(1,dimnull)).^[(L+1):(2*N)]*coefs;
    h = p - Nullmat*coefs; h= h';
else
    h = 1; 
end;
h = double(h);
if nargout==0
    disp('SG Filter impulse response (to be used with conv(x,h) command):');
    fliplr(h),
    disp('SG Filter expansion:');
    disp('      D^k ==> \nabla_k defined in the paper');
    disp('      S ==> convolution with all ones, i.e. S==>ones(1,2*N+1)');
end;
