function user_options = init_opt(self, control_mask, varargin)
% Initialize the optimization options and statistics.


fprintf('Optimization space dimension: %d\n', sum(sum(control_mask)));


%% statistics

self.stats.error = self.compute_error();
self.stats.wall_time = 0;
self.stats.cpu_time  = 0;
self.stats.integral  = self.seq.integral();


%% options

% MATLAB-style options processing.
% Converts a list of fieldname, value pairs in varargin to a struct.
user_options = struct(varargin{:});

% termination conditions and other options
defaults = struct(...
    'error_goal',        0.5 * (1e-3)^2 / self.system.norm2,...
    'max_loop_count',    1e10,...
    'max_walltime',      1800,...
    'max_cputime',       5e5,...
    'min_gradient_norm', 1e-20,...
    'plot_interval',     1);   % how often should we plot intermediate results?
[self.opt.options, user_options] = apply_options(defaults, user_options);


self.opt.initial_controls = self.seq.get();
self.opt.control_mask = control_mask;
self.opt.N_iter = 0;
self.opt.N_eval = 0;
self.opt.last_grad_norm = NaN;
self.opt.term_reason = 'None yet';
self.opt.stop = false;  % communication between the UI figure and monitor_func
self.opt.wall_start = now();
self.opt.cpu_start = cputime();
end


function [out, unused] = apply_options(defaults, opts)
% MATLAB-style options struct processing.
% Applies the options in the struct 'opts' to the struct 'defaults'.
% Returns the updated struct, and the struct of options that could not be parsed.

    % fields in opts
    names = fieldnames(opts);
    % logical array: are the corresponding fields present in defaults?
    present = isfield(defaults, names);
    
    % could not find a function for doing this
    out = defaults;
    for f = names(present).'
        out = setfield(out, f{1}, getfield(opts, f{1}));
    end
    % remove the fields we just used from opts
    unused = rmfield(opts, names(present));
end
