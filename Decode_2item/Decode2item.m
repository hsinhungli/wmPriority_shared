function [liks, hypers] = Decode2item(samples, p)

% TODO: ADD NOTES HERE

%% Initialize parameters
train_samples = samples(p.train_trials,:);
test_samples = samples(p.test_trials,:);
Ntesttrials = size(test_samples,1);
train_loc = p.stimpos(p.train_trials,:)/180*pi;
%test_ori = p.stimpos(p.test_trials,:)/180*pi;

high_loc = train_loc(:,1); %high priority stimulus is the same as the first column (target) in most trials
high_loc(p.condition(p.train_trials)==2) = train_loc(p.condition(p.train_trials)==2,2); %except the invalid cue trials
low_loc = train_loc(:,2);
low_loc(p.condition(p.train_trials)==2) = train_loc(p.condition(p.train_trials)==2,1);
clear samples

%%
Ctrain_1 = fun_basis(high_loc, p.nchan);
Ctrain_2 = fun_basis(low_loc, p.nchan);

n_angs = p.n_angs;
ang1 = p.ang1;
ang2 = p.ang2;
liks = nan(p.n_angs, p.n_angs, Ntesttrials);

%% Find best hyperparameter values (using inner CV-loop within the training data)
fprintf('\n--PERFORMING HYPERPARAMETER SEARCH\n')
lvr = linspace(0,1,50)';
lr = linspace(0,1,50)'; lr(1) = [];
hypers = find_lambda(p.runNs(p.train_trials), {lvr, lr});

%% estimate W & covariance
[W, est_w, noise] = estimate_W(train_samples, Ctrain_1, Ctrain_2, 0); %estimate_W(samples, C1, C2, do_boot, test_samples, test_C1, test_C2)
hypers(end+1) = est_w;

cov_est = estimate_cov(noise,hypers(1),hypers(2), W);

try
    prec_mat = invChol_mex(cov_est);
catch ME
    if strcmp(ME.identifier, 'MATLAB:invChol_mex:dpotrf:notposdef')
        fprintf('\nWARNING: Covariance estimate wasn''t positive definite in invChol_mex.\n');
    else
        rethrow(ME);
    end
    fprintf('Try use inv MATLAB function. \n');
    prec_mat = inv(cov_est);
end

%% Compute likelihood surface for test trials
for j = 1:Ntesttrials
    b = test_samples(j,:)';
    
    %[~, mll] = globminsearch; %Find the maximum likelihood so we can normalize the likelihood to a numerically feasible range
    % Compute the integral of the likelihood so we can normalize to a
    % probability distribution over s1,s2
    %Integ = integral2(@fun_lik,2*pi/1000,2*pi,2*pi/1000,2*pi);
    %lf(:,:,j) = reshape(fun_lik(ang1(:),ang2(:)),n_angs,n_angs)/Integ;
    
    this_lf = reshape(fun_lik(ang1(:),ang2(:)),n_angs,n_angs);
    liks(:,:,j) = this_lf/(sum(this_lf(:)));
end

fprintf \n

    function [minll, minder] = fun_negLL_norm(params)
        if nargout > 1
            [minll, minder] = fun_LL_norm(params);
        else
            minll = fun_LL_norm(params);
        end
        minll = -minll;
        if nargout>1
            minder = -minder;
        end
    end


    function [loglik, der] = fun_LL_norm(params)
        
        %Computes the log likelihood of the noise parameters. Also returns
        %the partial derivatives of each of the parameters (i.e. the
        %gradient), which are required by minimize.m and other efficient
        %optimization algorithms.
        
        nvox = size(noise,2);
        ntrials = size(noise,1);
        
        tau = params(1:end-2);
        sig = params(end-1);
        rho = params(end);
        
        
        [omi, NormConst] = invSNC(W(:,1:p.nchan), tau, sig, rho, p.singletau);
        
        XXt = noise'*noise;
        
        negloglik = 0.5*(MatProdTrace(XXt, omi) + ntrials*NormConst);
        
        if ~isreal(negloglik), negloglik = Inf; end %If we encounter a degenerate solution (indicated by complex-valued likelihoods), make sure that the likelihood goes to infinity.
        
        if any(vertcat(tau,sig)<0.001), negloglik = Inf; end
        if abs(rho) > 0.999999, negloglik = Inf; end
        
        loglik = -negloglik;
        
        if nargout > 1
            der = nan(size(params));
            
            ss = sqrt(ntrials);
            U = (omi*noise')/ss;
            
            dom = omi*(eye(nvox)-((1/ntrials)*XXt)*omi);
            
            JI = 1-eye(nvox);
            R = eye(nvox)*(1-rho) + rho;
            der(1:end-2) = 2*(dom.*R)*tau;
            der(end) = sum(sum(dom.*((tau*tau').*JI)));
            
            der(end-1) = 2*sig*MatProdTrace(W(:,1:p.nchan).'*omi, W(:,1:p.nchan)) - sum(sum((U.'*sqrt(2*sig)*W(:,1:p.nchan)).^2));
            
            der = -0.5*ntrials*der;
            
        end
        
        
    end

    function loss = fun_norm_loss(c_est, c0)
        
        try
            loss = (logdet(c_est, 'chol') + sum(sum(invChol_mex(c_est).*c0)))/size(c0,2);
        catch ME
            if any(strcmpi(ME.identifier, {'MATLAB:posdef', 'MATLAB:invChol_mex:dpotrf:notposdef'}))
                loss = (logdet(c_est) + trace(c_est\c0))/size(c0,2);
            else
                rethrow(ME);
            end
        end
        
    end

    function [W, est_w, noise, test_noise] = estimate_W(samples, C1, C2, do_boot, test_samples, test_C1, test_C2)
        if nargin <= 3, do_boot = 0; end
        N = size(C1,1);        
        
        if do_boot, idx = randi(N,N,1); else, idx = (1:N)'; end
        
        wlevel= p.wlevel;
        res = nan(length(wlevel),1);
        
        %Test out different gain factors for the low-priority item (just
        %find the minimum residual without considering cross-validation)
        for ww = 1:length(wlevel)
           thisC = C1 + wlevel(ww)*C2;
           thisW = thisC\samples;
           estB = thisC*thisW;
           res(ww) = sqrt(sum((samples(:) - estB(:)).^2));
        end
        
        %Can consider doing nested cv for this part, but didn't really
        %change the results...
        %ntrial = size(samples,1);
        %nfold = ntrial/12;
        %cvInd = kron((1:nfold)', ones(12,1));
        %for ww = 1:length(wlevel)
        %    res(ww) = 0;
        %    for cv_iter=1:nfold
        %        valInd = cvInd==cv_iter;
        %        thisC = C1 + wlevel(ww)*C2;
        %        thisW = thisC(~valInd,:)\samples(~valInd,:);
        %        estB = thisC(valInd,:)*thisW;
        %        res(ww) = res(ww) + sum(sum((samples(valInd,:) - estB).^2));
        %    end
        %end
        
        [~,Ind] = min(res);
        est_w = wlevel(Ind);
        estC = C1 + est_w*C2;
        W = (estC(idx,:)\samples(idx,:))';
        
        if nargout > 1
            noise = samples(idx,:) - estC(idx,:)*W';
            if nargout > 3
                test_C = test_C1 + est_w*test_C2;
                test_noise = test_samples - test_C*W';
            end
        end
    end

    function lambda = find_lambda(cvInd, lambda_range)
        % lambda_range can be a vector, in which case the same range is
        % used for the two lambda's. It can also be a matrix with two
        % columns, in which case the first column is used for lambda_var
        % and the second for lambda. Finally, it can be a cell array of
        % two cells, each of which contains a vector, in which case the
        % first vector is used for lambda_var, and the second for lambda
        % (this is the most flexible, as it allows you to specify two
        % different-length ranges for the two hyperparameters).
        
        if nargin<2, lambda_range = linspace(0,1,50)'; end
        
        cv_folds = unique(cvInd);
        K = length(cv_folds);
        
        assert(K>1, 'Must have at least two CV folds');
        
        clear W_cv est_noise_cv val_noise_cv
        W_cv{K} = []; est_noise_cv{K} = []; val_noise_cv{K} = [];
        
        % Pre-compute tuning weights and noise values to use in each
        % cross-validation split
        for cv_iter=1:K
            val_trials = cvInd==cv_folds(cv_iter);
            est_trials = ~val_trials;
            est_samples = train_samples(est_trials,:);
            val_samples = train_samples(val_trials,:);
            [W_cv{cv_iter}, ~, est_noise_cv{cv_iter}, val_noise_cv{cv_iter}] = estimate_W(est_samples, Ctrain_1(est_trials,:), Ctrain_2(est_trials,:), 0, val_samples, Ctrain_1(val_trials,:), Ctrain_2(val_trials,:));
        end
        
        % Grid search
        if ~iscell(lambda_range)
            if size(lambda_range,1)<size(lambda_range,2), lambda_range = lambda_range'; end
            if size(lambda_range,2)==1, lambda_range = repmat(lambda_range, 1,2); end
            lambda_range = mat2cell(lambda_range, size(lambda_range,1), [1 1]);
        end
        
        
        s = cellfun(@length, lambda_range);
        Ngrid = min(max(2, ceil(sqrt(s))), s); %Number of values to visit in each dimension (has to be at least 2, except if there is only 1 value for that dimension)
        
        grid_vec = cellfun(@(x,y) linspace(1, y, x), num2cell(Ngrid), num2cell(s), 'UniformOutput', 0);
        [grid_x, grid_y] = meshgrid(grid_vec{1}, grid_vec{2});
        [grid_l1, grid_l2] = meshgrid(lambda_range{1}, lambda_range{2});
        sz = fliplr(cellfun(@numel, lambda_range));
        
        fprintf('\n--GRID SEARCH--');
        losses = nan(numel(grid_x),1);
        for grid_iter=1:numel(grid_x)
            this_lambda = [lambda_range{1}(grid_x(grid_iter)) lambda_range{2}(grid_y(grid_iter))];
            losses(grid_iter) = visit(this_lambda);
            fprintf('\n %02d/%02d -- lambda_var: %3.2f, lambda: %3.2f, loss: %5.4g', [grid_iter, numel(grid_x), this_lambda, losses(grid_iter)]);
        end
        visited = sub2ind(sz, grid_y, grid_x); visited = visited(:);
        fprintf \n
        
        [best_loss, best_idx] = min(losses);
        best_idx = visited(best_idx);
        
        % Pattern search
        
        fprintf('\n--PATTERN SEARCH--');
        step_size = 2^floor(log2(diff(grid_y(1:2)/2))); %Round down to the nearest power of 2 (so we can keep dividing the step size in half)
        while 1
            [best_y,best_x] = ind2sub(sz, best_idx);
            new_x = best_x + [-1 1 -1 1]'*step_size;
            new_y = best_y + [-1 -1 1 1]'*step_size;
            del_idx = new_x<=0 | new_x> numel(lambda_range{1}) | new_y<=0 | new_y > numel(lambda_range{2});
            new_x(del_idx) = []; new_y(del_idx) = [];
            new_idx = sub2ind(sz, new_y, new_x);
            new_idx = new_idx(~ismember(new_idx, visited));
            if ~isempty(new_idx)
                this_losses = nan(size(new_idx));
                for ii = 1:length(new_idx)
                    this_lambda = [grid_l1(new_idx(ii)), grid_l2(new_idx(ii))];
                    this_losses(ii) = visit(this_lambda);
                    fprintf('\nStep size: %d, lambda_var: %3.2f, lambda: %3.2f, loss: %5.4g', [step_size, this_lambda, this_losses(ii)]);
                end
                visited = vertcat(visited, new_idx);
                losses = vertcat(losses, this_losses);
            end
            
            if any(this_losses<best_loss)
                [best_loss, best_idx] = min(losses);
                best_idx = visited(best_idx);
            elseif step_size>1
                step_size = step_size/2;
            else
                break
            end
        end
        fprintf \n
        
        lambda = [grid_l1(best_idx), grid_l2(best_idx)];
        
        fprintf('\nBest setting found: lambda_var = %3.2f, lambda = %3.2f, loss = %5.4g\n', [lambda, best_loss]);
        
        function loss = visit(lambda)
            loss = 0;
            for cv_iter2=1:K
                estC = estimate_cov(est_noise_cv{cv_iter2}, lambda(1), lambda(2), W_cv{cv_iter2});
                valC = (val_noise_cv{cv_iter2}'*val_noise_cv{cv_iter2})/size(val_noise_cv{cv_iter2},1); %sample covariance of validation data
                loss = loss + fun_norm_loss(estC, valC);
            end
            if imag(loss)~=0, loss = inf; end
        end
        
        
    end


    function C = estimate_cov(X,lambda_var,lambda, W)
        [n,pp] = size(X);
        W = W(:,1:p.nchan);
        
        vars = mean(X.^2);
        medVar = median(vars);
        
        t = tril(ones(pp),-1)==1;
        samp_cov = (X'*X/n);
        
        WWt = W*W';
        coeff = [WWt(t), ones(sum(t(:)),1)]\samp_cov(t); %lower left of the matrix and only off diagonal terms are included here
        
        target_diag = lambda_var*medVar + (1-lambda_var)*vars;
        target = coeff(1)*WWt + ones(pp)*coeff(2);
        target(eye(pp)==1)=target_diag;
        
        C = (1-lambda)*samp_cov + lambda*target;
        
        %ensure that the covariance matrix is well-conditioned
        [~, pp] = chol(C);
        orig_cond = cond(C);
        
        if pp>0 || orig_cond > 1e3
            [evec, eval] = eig(C);
            eval = diag(eval);
            max_eval = max(eval);
            eval = max(eval,max_eval/1000);
            C = evec*diag(eval)/evec;
            %fprintf('\nWARNING: Apply conditioning with max(condiiton number)=1000. \n');
        end
    end


    function out = MatProdTrace(mat1, mat2)
        % Computes the trace of a product of 2 matrices efficiently.
        mat2 = mat2';
        out = mat1(:)'*mat2(:);
    end

    function out = MatProdDiag(mat1, mat2)
        % Computes the diagonal of a product of 2 matrices efficiently.
        out = sum(mat1.*mat2', 2);
    end

    function out = fun_DKL(P, Q)
        % Computes KL-divergence from each row in Q to each corresponding
        % row in P
        z = P==0;
        out = P.*(log(P)-log(Q));
        out(z) = 0;
        out = sum(out,2);
    end

%% Functions used for 2D decoding
    function lik = fun_lik(s1,s2) %return likelihood
        %ll = exp(-fun_minLL([s1(:) s2(:)]) + mll);
        
        ll = -fun_minLL([s1(:) s2(:)]);
        lik = exp(ll - max(ll));
        lik = reshape(lik,size(s1,1),size(s1,2));
    end

    function negll = fun_minLL(s) %return negative log likelihood
        % W is nvox x nchan
        % b is nvox x 1
        % fun_basis(s)' is nchan x nangs
        
        s1 = s(:,1);
        s2 = s(:,2);
        
        bwc = b - W * (fun_basis(s1(:)) + est_w*fun_basis(s2(:)))'; %thisC = C1 + wlevel(ww)*C2;
        negll = 0.5*MatProdDiag(bwc'*prec_mat,bwc); 
    end

    function [sol, mll] = globminsearch
        % Finds the global minimum of the likelihood in s by doing a coarse
        % search in orientation space first.
        % (now does a 2-d search)
        [inits1,inits2] = meshgrid(linspace(0, 2*pi, 200),linspace(0, 2*pi, 200));
        inits1 = inits1(:); inits2 = inits2(:);
        fvals = fun_minLL([inits1 inits2]);
        [~, minI] = min(fvals);
        opts = optimset('MaxIter', 1e10, 'TolX', 1e-10, 'Display', 'off');
        [sol, mll] = fminsearch(@fun_minLL, [inits1(minI) inits2(minI)], opts);
    end


end
