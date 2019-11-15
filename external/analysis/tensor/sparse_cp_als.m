%%%%%%%%%%%%%%%%%%%%%%
%function to compute sparse cp-als
%requires Matlab Tensor Toolbox

% Reference:
%   Genevera I. Allen (2012). "Sparse Higher-Order Principal Components Analysis." 15th International Conference
%   on Artificial Intelligence and Statistics (AISTATS)

% Code retreived from author's website: http://www.stat.rice.edu/~gallen/software.html

function[U,V,W,D,bicmat,optlams] = sparse_cp_als(x,K,lamu,lamv,lamw,biclam,startu,startv,startw,posu,posv,posw,maxit)
ns = size(x);
n = ns(1); p = ns(2); q = ns(3);
r = length(lamu); setr = 1:r; setK = 1:K;
bicmat = []; optlams = [];
U = zeros(n,K,r); V = zeros(p,K,r); W = zeros(q,K,r); ds = zeros(K,r);
for j=1:r
    if norm(startu)==0
        U(:,:,j) = randn(n,K);
        V(:,:,j) = randn(p,K);
        W(:,:,j) = randn(q,K);
    else
        U(:,:,j) = startu;
        V(:,:,j) = startv;
        W(:,:,j) = startw;
    end    
    ind = 1; iter = 0; thr = 1e-6;
    while ind>thr & maxit>iter
        oldU = U(:,:,j); oldV = V(:,:,j); oldW = W(:,:,j);
        indi = 1; iteri = 0;
        krx = double(tenmat(x,1))*khatrirao(W(:,:,j),V(:,:,j));
        kri = krinnerprod(W(:,:,j),V(:,:,j));
        while indi>thr & iteri<maxit
            oldUi = U(:,:,j);
            for k=1:K
                U(:,k,j) = (1/kri(k,k))*soft_thr(krx(:,k) - U(:,setK~=k,j)*kri(setK~=k,k),lamu(j),posu);
            end
            iteri = iteri + 1; 
            indi = norm(U(:,:,j) - oldUi,'fro')/norm(oldUi,'fro');
        end
        ds(:,j) = max(abs(U(:,:,j)));
        if norm(ds(:,j))>0
            U(:,:,j) = U(:,:,j)*diag(1./ds(:,j));
        else
            U(:,:,j) = 0; V(:,:,j) = 0; W(:,:,j) = 0;
            break
        end
        
        
        indi = 1; iteri = 0;
        krx = double(tenmat(x,2))*khatrirao(W(:,:,j),U(:,:,j));
        kri = krinnerprod(W(:,:,j),U(:,:,j));
        while indi>thr & iteri<maxit
            oldVi = V(:,:,j);
            for k=1:K
                V(:,k,j) = (1/kri(k,k))*soft_thr(krx(:,k) - V(:,setK~=k,j)*kri(setK~=k,k),lamv(j),posv);
            end
            iteri = iteri + 1; 
            indi = norm(V(:,:,j) - oldVi,'fro')/norm(oldVi,'fro');
        end
        ds(:,j) = max(abs(V(:,:,j)));
        if norm(ds(:,j))>0
            V(:,:,j) = V(:,:,j)*diag(1./ds(:,j));
        else
            U(:,:,j) = 0; V(:,:,j) = 0; W(:,:,j) = 0;
            break
        end        

        indi = 1; iteri = 0;
        krx = double(tenmat(x,3))*khatrirao(V(:,:,j),U(:,:,j));
        kri = krinnerprod(V(:,:,j),U(:,:,j));
        while indi>thr & iteri<maxit
            oldWi = W(:,:,j);
            for k=1:K
                W(:,k,j) = (1/kri(k,k))*soft_thr(krx(:,k) - W(:,setK~=k,j)*kri(setK~=k,k),lamw(j),posw);
            end
            iteri = iteri + 1; 
            indi = norm(W(:,:,j) - oldWi,'fro')/norm(oldWi,'fro');
        end
        ds(:,j) = max(abs(W(:,:,j)));
        if norm(ds(:,j))>0
            W(:,:,j) = W(:,:,j)*diag(1./ds(:,j));
        else
            U(:,:,j) = 0; V(:,:,j) = 0; W(:,:,j) = 0;
            break
        end        

        ind = norm(oldU - U(:,:,j),'fro')/norm(oldU,'fro') +norm(oldV - V(:,:,j),'fro')/norm(oldV,'fro')+norm(oldW - W(:,:,j),'fro')/norm(oldW,'fro'); 
        iter = iter + 1;                
    end    
    %calculate correct Ds
    %order U's according to D's
    for i=1:K
        U(:,i,j) = U(:,i,j)/norm(U(:,i,j));
        V(:,i,j) = V(:,i,j)/norm(V(:,i,j));
        W(:,i,j) = W(:,i,j)/norm(W(:,i,j));
    end
    Dtmp = ttm(ttm(ttm(x,U(:,:,j)',1),V(:,:,j)',2),W(:,:,j)',3);
    for i=1:K
        Ds(i) = Dtmp(i,i,i);
    end    
    [dtmp,dind] = sort(Ds,'descend');
    ds(:,j) = dtmp;
    U(:,:,j) = U(:,dind,j);
    V(:,:,j) = V(:,dind,j);
    W(:,:,j) = W(:,dind,j);
    switch biclam
      case 'u'
        df = sum(sum(U(:,:,j)~=0)); 
      case 'v'
        df = sum(sum(V(:,:,j)~=0));
      case 'w'
        df = sum(sum(W(:,:,j)~=0));
    end
    xr = sum(sum(sum(double((x - full(ktensor(ds(:,j),U(:,:,j),V(:,:,j),W(:,:,j)))).^2))));
    bicmat(j) = log( xr / (n*p*q) ) + (log(n*p*q)/(n*p*q)).*df;                
end
if r==1
    U = U(:,:,1); V = V(:,:,1); W = W(:,:,1); D = ds(:,1);
else
    ind = bicmat==min(bicmat);
    if sum(ind)>1
        ind = min(setr(ind));
    end
    switch biclam
      case 'u'
    optlams = [optlams; lamu(ind)];
      case 'v'
        optlams = [optlams; lamv(ind)];
      case 'w'
        optlams = [optlams; lamw(ind)];
    end
    U = U(:,:,ind); V = V(:,:,ind); W = W(:,:,ind);
    D = ds(:,ind);
end



