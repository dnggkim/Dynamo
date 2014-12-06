function term_reason = search_BFGS(self, control_mask, user_options)
% BFGS optimization.


if nargin < 3
    user_options = {};
end

% initialize the optimization data structures
self.init_opt(control_mask, user_options);

% default options for fminunc
problem.options = optimset(...
    'MaxIter',      1e4,...
    'TolX',         1e-8,...
    'TolFun',       1e-8,...
    'DerivativeCheck', 'off',...
    'FinDiffType',  'central',...    
    'GradObj',      'on',...
    'LargeScale',   'off', ...
    'OutputFcn',    @(x, optimValues, state) monitor_func(self, x, optimValues, state),...
    'Display',      'off');

% additional user-defined options
problem.options = optimset(problem.options, user_options);

% save a copy
self.opt.BFGS_options = problem.options;

% define the optimization problem
problem.objective = @(x) goal_and_gradient_function_wrapper(self, x);
problem.x0 = self.seq.get(self.opt.control_mask);
problem.solver = 'fminunc';

fprintf('\nOptimizing algorithm: BFGS. Running...\n\n'); drawnow;

% try to minimise objective function to zero
[x, cost, exitflag, output] = fminunc(problem.objective, problem.x0, problem.options);

self.update_controls(x, self.opt.control_mask); % It may be different than the last point evaluated
term_reason = self.opt.term_reason;
end


function [err, grad] = goal_and_gradient_function_wrapper(self, x)
% x is a vector containing (a subset of) the controls
    
    self.opt.N_eval = self.opt.N_eval + 1;

    self.update_controls(x, self.opt.control_mask);
    [err, grad] = self.compute_error(self.opt.control_mask);
    self.opt.last_grad_norm = sqrt(sum(sum(grad .* grad)));
end
