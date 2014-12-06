function init_opt(self, control_mask, options)
% Initialize the optimization options and statistics.


%% options

self.opt.initial_controls = self.seq.get();
self.opt.control_mask = control_mask;


fprintf('Optimization space dimension: %d\n', sum(sum(control_mask)));

self.opt.N_iter = 0;
self.opt.N_eval = 0;
self.opt.last_grad_norm = NaN;
self.opt.term_reason = 'none yet';


error_goal = 0.5 * (1e-4)^2 / self.system.norm2


self.opt.term_cond = struct( ...
    'max_loop_count',     1e10, ...
    'error_goal',         error_goal, ...
    'max_wall_time',       1800, ...
    'max_cputime',         5e5, ...
    'min_gradient_norm', 1e-20);

if isfield(options, 'max_wall_time')
    % TEST FIXME NOW
    self.opt.term_cond.max_wall_time = options.max_wall_time
end


% communication between the UI figure and monitor_func
self.opt.stop = false;

% should we plot intermediate results?
if isfield(options, 'plot_interval')
    self.opt.plot_interval = options.plot_interval;
else
    self.opt.plot_interval = 0;
end

self.opt.wall_start = now();
self.opt.cpu_start = cputime();


%% statistics
self.stats.error = self.compute_error();
self.stats.wall_time = 0;
self.stats.cpu_time  = 0;
self.stats.integral  = self.seq.integral();
end
