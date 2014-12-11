function term_reason = search_BFGS(self, control_mask, varargin)
% BFGS optimization.


if nargin < 2
    % by default update all the controls except taus
    control_mask = self.full_mask(false);
end

% common initialization of the optimization data structures
user_options = self.init_opt(control_mask, varargin{:});

% define the optimization problem
problem.objective = @(x) goal_and_gradient_function_wrapper(self, x);
problem.x0 = self.seq.get(self.opt.control_mask);
problem.solver = 'fminunc';

% default options for fminunc

% TODO with newer MATLAB versions we would do it like this:
%problem.options = optimoptions('fminunc',...
%    'Algorithm',    'quasi-newton',...
% for now, use the old optimset()
problem.options = optimset(...
    'DerivativeCheck', 'off',...
    'Display',         'final',...
    'GradObj',         'on',...    % use user-supplied gradient
    'LargeScale',      'off', ...  % force quasi-newton algorithm (BFGS)
    'MaxIter',         1e4,...
    'OutputFcn', @(x, optimValues, state) monitor_func(self, x, optimValues, state),...
    'TolFun',          1e-8,...
    'TolX',            1e-8);
% additional user-defined options
problem.options = optimset(problem.options, user_options);
% save a copy
self.opt.BFGS_options = problem.options;

fprintf('\nOptimizing algorithm: BFGS. Running...\n\n'); drawnow;

% try to minimise objective function to zero
[x, cost, exitflag, output] = fminunc(problem);

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
