
% Example call for tensor X with rank R:
% cprand(X, R, 'num_samples', 50, 'printitn', 50, 'init', 'random', 'maxiters', MAX_ITERS, 'desired_fit', TOL, 'fft', 1); 

function [P,Uinit,output] = cprand(X,R,varargin)

    N = ndims(X);
    normX = norm(X);
    sz = size(X);

    %% Set algorithm parameters from input or by using defaults
    params = inputParser;
    params.addParamValue('tol',1e-4,@isscalar); %used to be 1e-5
    params.addParamValue('maxiters',100,@(x) isscalar(x) & x > 0);
    params.addParamValue('dimorder',1:N,@(x) isequal(sort(x),1:N));
    params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
    params.addParamValue('printitn',10,@isscalar);
    params.addParamValue('desired_fit',0.99,@(x) isscalar(x) & x > 0 & x < 1);
    params.addParamValue('fft',0,@(x) isscalar(x) & x == 0 || x == 1);
    params.addParamValue('num_samples',0,@(x) isscalar(x) & x > 0);
    params.addParamValue('stopping',1,@(x) isscalar(x) & x == 0 || x == 1);

    params.parse(varargin{:});

    %% Copy from params object
    fitchangetol = params.Results.tol;
    maxiters = params.Results.maxiters;
    dimorder = params.Results.dimorder;
    init = params.Results.init;
    printitn = params.Results.printitn;
    desired_fit = params.Results.desired_fit;
    do_fft = params.Results.fft;
    stopping = params.Results.stopping;

    % Set number of samples.
    % Setting samples too low may result in rank-deficiency and a warning will appear.
    if (params.Results.num_samples > 0)
        num_samples = params.Results.num_samples;
    else
        num_samples = ceil(10*R*log2(R));
        fprintf('Warning: num_samples not passed to CPRAND.\n Using (10 * R logR)=%d as default.\n',num_samples); 
    end

    %% Error checking 
    %% Set up and error checking on initial guess for U.
    if iscell(init)
        Uinit = init;
        if numel(Uinit) ~= N
            error('OPTS.init does not have %d cells',N);
        end
        for n = dimorder(2:end);
            if ~isequal(size(Uinit{n}),[size(X,n) R])
                error('OPTS.init{%d} is the wrong size',n);
            end
        end
    else
        % Observe that we don't need to calculate an initial guess for the
        % first index in dimorder because that will be solved for in the first
        % inner iteration.
        if strcmp(init,'random')
            Uinit = cell(N,1);
            for n = dimorder(2:end)
                Uinit{n} = rand(sz(n),R);
            end
        elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
            Uinit = cell(N,1);
            for n = dimorder(2:end)
                Uinit{n} = nvecs(X,n,R);
            end
        else
            error('The selected initialization method is not supported');
        end
    end

    %% Set up for iterations - initializing U and the fit.
    U = Uinit;
    fit = 0;
    diag_flips = [];

    if printitn>0
        if(do_fft)
            fprintf('CP-RAND-FFT: ');
        else
            fprintf('CP_RAND: ');
        end
    end

    %% Sample input tensor for stopping criterion
    fitsamples = min(nnz(X),2^14);
    [Xfit_subs, Xfit_idxs] = sample_all_modes(fitsamples, sz);
    Xfit_vals = X(Xfit_subs);
    Xfit = sptensor(Xfit_subs,Xfit_vals,sz);
    normXf = norm(Xfit);
    %% Main Loop: Iterate until convergence
    if (do_fft)
        % Compute random diagonal D_n for each factor
        diag_flips = cell(N,1);
        for n = 1:N
            diag_flips{n} = (rand(sz(n),1)<0.5)*2-1;
        end

        X_mixed = X.data;
        % Mixing is equivalent to a series of TTMs with D_n, F_n
        % However, we can use bsxfun and fft to avoid matricizing.
        for n = N:-1:1
            % Reshape the diagonal flips into a 1*...*sz(n)*...*1 tensor
            % This lets us use bsxfun along the nth dimension.
            bsxdims = ones(1,N);
            bsxdims(n) = sz(n);
            flips = reshape(diag_flips{n},bsxdims);
            % fft(...,[],n) operates fiber-wise on dimension n
            X_mixed = fft(bsxfun(@times, X_mixed, flips),[],n);
        end
        X_mixed = tensor(X_mixed);
    else
        X_mixed = X; % no mixing
    end

    % ALS Loop
    for iter = 1:maxiters
        fitold = fit;
        % Iterate over all N modes of the tensor
        for n = dimorder(1:end)

            mixinfo.dofft = do_fft;
            mixinfo.signflips = diag_flips;
            [Unew, sampX, sampZ]= dense_sample_mttkrp(X_mixed,U,n,num_samples,mixinfo);

            if issparse(Unew)
                Unew = full(Unew);   % for the case R=1
            end

            % Normalize each vector to prevent singularities in coefmatrix
            if iter == 1
                lambda = sqrt(sum(Unew.^2,1))'; %2-norm
            else
                lambda = max( max(abs(Unew),[],1), 1 )'; %max-norm
            end      

            Unew = bsxfun(@rdivide, Unew, lambda');
            U{n} = Unew;
        end

        P = ktensor(lambda, U);
        if (stopping || mod(iter,printitn)==0)
            if normX == 0
                fit = norm(P)^2 - 2 * innerprod(X,P);
            else
                Ps = sample_ktensor(P, Xfit_subs);
                diff_mean = mean((Xfit_vals - Ps).^2);
                fit = 1 - sqrt(diff_mean*prod(size(X)))/normX;
            end
            fitchange = abs(fitold - fit);
            %SCP convergence criteria
            if (stopping) && (iter > 1) && ((fitchange < fitchangetol) || fit > desired_fit)
                flag = 0;
            else
                flag = 1;
            end
            
            if (mod(iter,printitn)==0) %|| ((printitn>0) && (flag==0))
                fprintf(' Iter %2d: f = %e f-delta = %7.1e\n', iter, fit, fitchange); 
            end
            
            % Check for convergence
            if (flag == 0)
                break;
            end        
        end
    end   
    %% Clean up final result
    % Arrange the final tensor so that the columns are normalized.
    P = arrange(P);
    P = fixsigns(P); % Fix the signs

    %if printitn>0
    if normX == 0
        fit = norm(P)^2 - 2 * innerprod(X,P);
    else
        normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
        fit = 1 - (normresidual / normX);%fraction explained by model
        Ps = sample_ktensor(P, Xfit_subs);
        Xfit_mean = mean((Xfit_vals - Ps).^2);
        testfit = 1 - sqrt(Xfit_mean*prod(size(X)))/normX;
    end
    if printitn>0
        fprintf(' Final fit = %e Final estimated fit = %e \n', fit, testfit);
    end
    %end

    output = struct;
    output.params = params.Results;
    output.iters = iter;
    output.fit = fit; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MTTKRP That performs sampling after transforming X and the KR-product with an FFT
% And then solves using normal equations
function [V, Xsamp, Zsamp] = dense_sample_mttkrp(X,U,n,num_samples,mixinfo)
    N = ndims(X);
    if (N < 2) error('MTTKRP is invalid for tensors with fewer than 2 dimensions'); end
    if (length(U) ~= N) error('Cell array is the wrong length'); end

    if n == 1 
        R = size(U{2},2);
    else 
        R = size(U{1},2); 
    end

    for i = 1:N
        if i == n, continue; end
        if (size(U{i},1) ~= size(X,i)) || (size(U{i},2) ~= R)
            error('Entry %d of cell array is wrong size', i);
        end
    end
  
    dims = size(X);

    if (mixinfo.dofft)
        % Mix factor matrices: U{i} = F{i}*D{i}*U{i}
        for i = [1:n-1,n+1:N]
            U{i} = fft(bsxfun(@times,U{i},mixinfo.signflips{i}));
        end
    end

    % Compute uniform samples for tensor and factor matrices
    [tensor_idx, factor_idx] = sample_mode_n(num_samples, dims, n);

    % Reshape the sampled tensor
    Xsamp = reshape(X(tensor_idx), dims(n), []);

    % Perform a sampled KRP
    Zsamp = skr(U{[1:n-1,n+1:N]}, factor_idx);

    if (mixinfo.dofft)
        % Preprocess the sampled tensor: Xsamp = D^-1 F^-1 Xsamp
        Xsamp = bsxfun(@times, ifft(Xsamp), mixinfo.signflips{n}); 
    end

    %Solve step
    XX = [real(Xsamp) imag(Xsamp)];
    ZZ = [real(Zsamp); imag(Zsamp)];
    V = XX / ZZ.';
    % V = real(Xsamp / Zsamp.');
    return;
end


% Sample Khatri-Rao Product of a cell array of factors
% Without forming the full KR Product
function P = skr(varargin)
    if iscell(varargin{1}) % Input is a single cell array
        A = varargin{1};
    else % Input is a sequence of matrices
        A = {varargin{1:end-1}};
    end

    numfactors = length(A);
    matorder = numfactors:-1:1;
    idxs = varargin{end};

    %% Error check on matrices and compute number of rows in result 
    ndimsA = cellfun(@ndims, A);
    if(~all(ndimsA == 2))
        error('Each argument must be a matrix');
    end

    ncols = cellfun(@(x) size(x, 2), A);
    if(~all(ncols == ncols(1)))
        error('All matrices must have the same number of columns.');
    end

    P = A{matorder(1)}(idxs(:,matorder(1)),:);
    for i = matorder(2:end)
        %P = P .*A{i}(idxs(:,i),:);
        P = bsxfun(@times, P, A{i}(idxs(:,i),:));
    end
end


% Random sample fibers in mode n from tensor X
% Generate the corresponding indices for the factor matrices as a tuple
function [tensor_idx, factor_idx] = sample_mode_n(num_samples, dims, n)
    D = length(dims);
    tensor_idx = zeros(num_samples, D);     % Tuples that index fibers in original tensor

    tensor_idx(:,n) = ones(num_samples, 1);
    for i = [1:n-1,n+1:D]
        % Uniformly sample w.r. in each dimension besides n
        tensor_idx(:,i) = randi(dims(i), num_samples, 1);
    end

    % Save indices to sample from factor matrices
    factor_idx = tensor_idx(:,[1:n-1,n+1:D]);

    % Expand tensor_idx so that every fiber element is included
    %tensor_idx = repelem(tensor_idx,dims(n),1); % not portable
    tensor_idx = kron(tensor_idx,ones(dims(n),1)); % portable
    tensor_idx(:,n) = repmat([1:dims(n)]',num_samples,1);
    tensor_idx = tt_sub2ind(dims, tensor_idx);
end


% Random sample fibers in mode n from tensor X
% Generate the corresponding indices for the factor matrices as a tuple
function [subs, idxs] = sample_all_modes(num_samples, dims)
    D = length(dims);
    subs = zeros(num_samples, D);     % Tuples that index fibers in original tensor

    for i = [1:D]
        % Uniformly sample w.r. in each dimension
        subs(:,i) = randi(dims(i), num_samples, 1);
    end

    % subs can be used to sample from factor matrices as well as the tensor
    subs = unique(subs,'rows'); %todo: do this more efficiently
    idxs = tt_sub2ind(dims, subs);
end


% Random sample fibers in mode n from tensor X
% Generate the corresponding indices for the factor matrices as a tuple
function [data] = sample_ktensor(P, subs)
    sz = size(P);
    data = skr(P.u, subs) * P.lambda;
    %Xk = sptensor(subs,data,sz);
    % we may rarely see a 0 value in data
    % so we should just return data
end






% % Example call for tensor X with rank R:
% % cprand(X, R, 'num_samples', 50, 'printitn', 50, 'init', 'random', 'maxiters', MAX_ITERS, 'desired_fit', TOL, 'fft', 1); 

% function [P,Uinit,output,X_mixed,Q] = cprand(X,R,varargin)

%     N = ndims(X);
%     normX = norm(X);
%     sz = size(X);

%     %% Set algorithm parameters from input or by using defaults
%     params = inputParser;
%     params.addParamValue('tol',1e-5,@isscalar); %used to be 1e-5
%     params.addParamValue('maxiters',100,@(x) isscalar(x) & x > 0);
%     params.addParamValue('dimorder',1:N,@(x) isequal(sort(x),1:N));
%     params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
%     params.addParamValue('printitn',10,@isscalar);
%     params.addParamValue('desired_fit',0.99,@(x) isscalar(x) & x > 0 & x < 1);
%     params.addParamValue('fft',0,@(x) isscalar(x) & x == 0 || x == 1);
%     params.addParamValue('num_samples',0,@(x) isscalar(x) & x > 0);
%     params.addParamValue('stopping',1,@(x) isscalar(x) & x == 0 || x == 1);

%     params.parse(varargin{:});

%     %% Copy from params object
%     fitchangetol = params.Results.tol;
%     maxiters = params.Results.maxiters;
%     dimorder = params.Results.dimorder;
%     init = params.Results.init;
%     printitn = params.Results.printitn;
%     desired_fit = params.Results.desired_fit;
%     do_fft = params.Results.fft;
%     stopping = params.Results.stopping;

%     % Set number of samples.
%     % Setting samples too low may result in rank-deficiency and a warning will appear.
%     if (params.Results.num_samples > 0)
%         num_samples = params.Results.num_samples;
%     else
%         num_samples = ceil(10*R*log2(R));
%         fprintf('Warning: num_samples not passed to CPRAND.\n Using (10 * R logR)=%d as default.\n',num_samples); 
%     end

%     diag_flips = [];
%     if (do_fft)
%         % Compute random diagonal D_n for each factor
%         diag_flips = cell(N,1);
%         for n = 1:N
%             diag_flips{n} = (rand(sz(n),1)<0.5)*2-1;
%         end
%     end

%     %% Error checking 
%     %% Set up and error checking on initial guess for U.
%     if iscell(init)
%         Uinit = init;
%         if numel(Uinit) ~= N
%             error('OPTS.init does not have %d cells',N);
%         end
%         for n = dimorder(2:end);
%             if ~isequal(size(Uinit{n}),[size(X,n) R])
%                 error('OPTS.init{%d} is the wrong size',n);
%             end
%         end
%     else
%         % Observe that we don't need to calculate an initial guess for the
%         % first index in dimorder because that will be solved for in the first
%         % inner iteration.
%         if strcmp(init,'random')
%             Uinit = cell(N,1);
%             for n = dimorder(2:end)
%                 Uinit{n} = rand(sz(n),R);
%                 if (do_fft)
%                     Uinit{n} = fft(bsxfun(@times,Uinit{n},diag_flips{n}));
%                 end
%             end
%         elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
%             Uinit = cell(N,1);
%             for n = dimorder(2:end)
%                 Uinit{n} = nvecs(X,n,R);
%                 if (do_fft)
%                     Uinit{n} = fft(bsxfun(@times,Uinit{n},diag_flips{n}));
%                 end
%             end
%         else
%             error('The selected initialization method is not supported');
%         end
%     end

%     %% Set up for iterations - initializing U and the fit.
%     U = Uinit;
%     fit = 0;

%     if printitn>0
%         if(do_fft)
%             fprintf('CP-RAND-FFT: ');
%         else
%             fprintf('CP_RAND: ');
%         end
%     end

%     %% Main Loop: Iterate until convergence
%     if (do_fft)
%         X_mixed = X.data;
%         % Mixing is equivalent to a series of TTMs with D_n, F_n
%         % However, we can use bsxfun and fft to avoid matricizing.
%         for n = N:-1:1
%             % Reshape the diagonal flips into a 1*...*sz(n)*...*1 tensor
%             % This lets us use bsxfun along the nth dimension.
%             bsxdims = ones(1,N);
%             bsxdims(n) = sz(n);
%             flips = reshape(diag_flips{n},bsxdims);

%             % fft(...,[],n) operates fiber-wise on dimension n
%             X_mixed = fft(bsxfun(@times, X_mixed, flips),[],n);
%         end
%         X_mixed = tensor(X_mixed);
%     else
%         X_mixed = X; % no mixing
%     end
%     normXm = norm(X_mixed);
% %if (do_fft)
%     %% Sample input tensor for stopping criterion
%     fitsamples = min(nnz(X_mixed),1024);
%     [Xfit_subs, Xfit_idxs] = sample_all_modes(fitsamples, sz);
%     Xfit_vals = X(Xfit_subs);
%     Xfit = sptensor(Xfit_subs,Xfit_vals,sz);
%     normXf = norm(Xfit);
%         %normX = norm(X_mixed);
%     % else
%     %     %% Sample input tensor for stopping criterion
%     %     fitsamples = min(nnz(X),1024);
%     %     [Xfit_subs, Xfit_idxs] = sample_all_modes(fitsamples, sz);
%     %     Xfit_vals = X_mixed(Xfit_subs);
%     %     Xfit = sptensor(Xfit_subs,Xfit_vals,sz);
%     %     normXf = norm(Xfit);
%     %     normX = norm(X_mixed);
%     % end

%     % ALS Loop
%     for iter = 1:maxiters
%         fitold = fit;
%         % Iterate over all N modes of the tensor
%         for n = dimorder(1:end)
%             mixinfo.dofft = do_fft;
%             mixinfo.signflips = diag_flips;
%             [Unew, sampX, sampZ]= dense_sample_mttkrp(X_mixed,U,n,num_samples,mixinfo);

%             if issparse(Unew)
%                 Unew = full(Unew);   % for the case R=1
%             end

%             %Normalize each vector to prevent singularities in coefmatrix
%             if iter == 1
%                 lambda = sqrt(sum(Unew.^2,1))'; %2-norm
%             else
%                 lambda = max( max(abs(Unew),[],1), 1 )'; %max-norm
%             end
%             Unew = bsxfun(@rdivide, Unew, lambda');
%             U{n} = Unew;
%         end
%         P = ktensor(lambda, U);
%         %P = ktensor(U);
%         if (stopping || mod(iter,printitn)==0)
%             if do_fft
%                 if normX == 0
%                     fit = norm(P)^2 - 2 * innerprod(X_mixed,P);
%                 else
%                     Q = P;
%                     for i = [1:N]
%                         Q{i} = ifft(Q{i});
%                         Q{i} = bsxfun(@times,Q{i},diag_flips{i});
%                     end
%                     Q = arrange(Q);
%                     Q = fixsigns(Q); 
%                     normresidual = sqrt( normX^2 + norm(Q)^2 - 2 * innerprod(X,Q) );
%                     fit = 1 - (normresidual / normX);%fraction explained by model
%                     %fit = norm(X_mixed - full(Q)) / norm(X_mixed)
%                     % normresidual = sqrt( normXm^2 + norm(P)^2 - 2 * innerprod(X_mixed,P) );
%                     % fit = 1 - (normresidual / normXm);%
%                     % Ps = sample_ktensor(P, Xfit_subs);
%                     % diff_mean = mean((Xfit_vals - Ps).^2)
%                     % fit = 1 - sqrt(diff_mean*prod(size(X_mixed)))/normX
%                 end
%             else
%                 if normX == 0
%                     fit = norm(P)^2 - 2 * innerprod(X,P);
%                 else
%                     Ps = sample_ktensor(P, Xfit_subs);
%                     diff_mean = mean((Xfit_vals - Ps).^2);
%                     fit = 1 - sqrt(diff_mean*prod(size(X)))/normX;
%                 end 
%             end
%             fitchange = abs(fitold - fit);
%             % SCP convergence criteria
%             if (stopping) && (iter > 1) && ((fitchange < fitchangetol) || fit > desired_fit)
%                 flag = 0;
%             else
%                 flag = 1;
%             end
            
%             if (mod(iter,printitn)==0) %|| ((printitn>0) && (flag==0))
%                 fprintf(' Iter %2d: f = %e f-delta = %7.1e\n', iter, fit, fitchange); 
%             end
            
%             % Check for convergence
%             if (flag == 0)
%                 break;
%             end        
%         end
%     end   
%     %% Clean up final result
%     % Arrange the final tensor so that the columns are normalized.
%     if (do_fft)
%         for i = [1:N]
%             P{i} = ifft(P{i});
%             P{i} = bsxfun(@times,P{i},diag_flips{i});
%         end
%     end
%     Q = P;
%     P = arrange(P);
%     P = fixsigns(P); % Fix the signs
%     %if printitn>0
%     if normX == 0
%         fit = norm(P)^2 - 2 * innerprod(X,P);
%     else
%         %normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
%         %fit = 1 - (normresidual / normX);%fraction explained by model
%         fit = norm(X - full(P))/norm(X);
%         Ps = sample_ktensor(P, Xfit_subs);
%         Xfit_mean = mean((Xfit_vals - Ps).^2);
%         testfit = 1 - sqrt(Xfit_mean*prod(size(X)))/normX;
%     end
%     fprintf(' Final fit = %e Final estimated fit = %e \n', fit, testfit);
%     %end

%     output = struct;
%     output.params = params.Results;
%     output.iters = iter;
%     output.fit = fit; 
% end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%% Sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % MTTKRP That performs sampling after transforming X and the KR-product with an FFT
% % And then solves using normal equations
% function [V, Xsamp, Zsamp] = dense_sample_mttkrp(X,U,n,num_samples,mixinfo)
%     N = ndims(X);
%     if (N < 2) error('MTTKRP is invalid for tensors with fewer than 2 dimensions'); end
%     if (length(U) ~= N) error('Cell array is the wrong length'); end

%     if n == 1 
%         R = size(U{2},2);
%     else 
%         R = size(U{1},2); 
%     end

%     for i = 1:N
%         if i == n, continue; end
%         if (size(U{i},1) ~= size(X,i)) || (size(U{i},2) ~= R)
%             error('Entry %d of cell array is wrong size', i);
%         end
%     end
  
%     dims = size(X);

%     % if (mixinfo.dofft)
%     %     % Mix factor matrices: U{i} = F{i}*D{i}*U{i}
%     %     for i = [1:n-1,n+1:N]
%     %         U{i} = fft(bsxfun(@times,U{i},mixinfo.signflips{i}));
%     %     end
%     % end

%     % Compute uniform samples for tensor and factor matrices
%     [tensor_idx, factor_idx] = sample_mode_n(num_samples, dims, n);

%     % Reshape the sampled tensor
%     Xsamp = reshape(X(tensor_idx), dims(n), []);

%     % Perform a sampled KRP
%     Zsamp = skr(U{[1:n-1,n+1:N]}, factor_idx);

%     % if (mixinfo.dofft)
%     %     % Preprocess the sampled tensor: Xsamp = D^-1 F^-1 Xsamp
%     %     Xsamp = bsxfun(@times, ifft(Xsamp), mixinfo.signflips{n}); 
%     % end

%     % Solve step
%     V = Xsamp / Zsamp.';
%     return;
% end


% % Sample Khatri-Rao Product of a cell array of factors
% % Without forming the full KR Product
% function P = skr(varargin)
%     if iscell(varargin{1}) % Input is a single cell array
%         A = varargin{1};
%     else % Input is a sequence of matrices
%         A = {varargin{1:end-1}};
%     end

%     numfactors = length(A);
%     matorder = numfactors:-1:1;
%     idxs = varargin{end};

%     %% Error check on matrices and compute number of rows in result 
%     ndimsA = cellfun(@ndims, A);
%     if(~all(ndimsA == 2))
%         error('Each argument must be a matrix');
%     end

%     ncols = cellfun(@(x) size(x, 2), A);
%     if(~all(ncols == ncols(1)))
%         error('All matrices must have the same number of columns.');
%     end

%     P = A{matorder(1)}(idxs(:,matorder(1)),:);
%     for i = matorder(2:end)
%         %P = P .*A{i}(idxs(:,i),:);
%         P = bsxfun(@times, P, A{i}(idxs(:,i),:));
%     end
% end


% % Random sample fibers in mode n from tensor X
% % Generate the corresponding indices for the factor matrices as a tuple
% function [tensor_idx, factor_idx] = sample_mode_n(num_samples, dims, n)
%     D = length(dims);
%     tensor_idx = zeros(num_samples, D);     % Tuples that index fibers in original tensor

%     tensor_idx(:,n) = ones(num_samples, 1);
%     for i = [1:n-1,n+1:D]
%         % Uniformly sample w.r. in each dimension besides n
%         tensor_idx(:,i) = randi(dims(i), num_samples, 1);
%     end

%     % Save indices to sample from factor matrices
%     factor_idx = tensor_idx(:,[1:n-1,n+1:D]);

%     % Expand tensor_idx so that every fiber element is included
%     %tensor_idx = repelem(tensor_idx,dims(n),1); % not portable
%     tensor_idx = kron(tensor_idx,ones(dims(n),1)); % portable
%     tensor_idx(:,n) = repmat([1:dims(n)]',num_samples,1);
%     tensor_idx = tt_sub2ind(dims, tensor_idx);
% end


% % Random sample fibers in mode n from tensor X
% % Generate the corresponding indices for the factor matrices as a tuple
% function [subs, idxs] = sample_all_modes(num_samples, dims)
%     D = length(dims);
%     subs = zeros(num_samples, D);     % Tuples that index fibers in original tensor

%     for i = [1:D]
%         % Uniformly sample w.r. in each dimension
%         subs(:,i) = randi(dims(i), num_samples, 1);
%     end

%     % subs can be used to sample from factor matrices as well as the tensor
%     subs = unique(subs,'rows'); %todo: do this more efficiently
%     idxs = tt_sub2ind(dims, subs);
% end


% % Random sample fibers in mode n from tensor X
% % Generate the corresponding indices for the factor matrices as a tuple
% function [data] = sample_ktensor(P, subs)
%     sz = size(P);
%     data = skr(P.u, subs) * P.lambda;
%     %Xk = sptensor(subs,data,sz);
%     % we may rarely see a 0 value in data
%     % so we should just return data
% end



