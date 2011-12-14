function dynamo_init_opt(control_mask)
% Initialize the optimization module.
    
global OC;

OC.opt.control_mask = control_mask;

fprintf('Optimization space dimension: %d\n', sum(sum(control_mask)));

OC.opt.N_iter = 0;
OC.opt.N_eval = 0;
OC.opt.last_grad_norm = NaN;
OC.opt.term_reason = 'none yet';


OC.opt.term_cond = struct( ...
    'max_loop_count',     1e10, ...
    'error_goal',        1e-10, ...
    'max_wall_time',       180, ...
    'max_cputime',         180, ...
    'min_gradient_norm', 1e-20);


OC.opt.wall_start = now();
OC.opt.cpu_start = cputime();


OC.stats.error = [];
OC.stats.wall_time = [];
OC.stats.cpu_time  = [];
OC.stats.fluence  = [];
