function analyze(self)
% Analyzes the results of an optimization run.


err = self.compute_error();
fprintf('Final normalized error: %g\n    Wall time: %g s\n    CPU  time: %g s\nTermination reason: %s\n\n\n', ...
	err, self.stats.wall_time(end), self.stats.cpu_time(end), self.opt.term_reason);

fprintf('Number of gradient evaluations: %d\n', self.opt.N_eval);
fprintf('Final sequence duration: %g\n', sum(self.seq.tau));


% rough error estimates
%d = real(eig(self.system.A));
%n = length(d);
%T = sum(self.seq.tau);
%e_max = 1-exp(-T*sum(d)/n)
%e_min = 1-sum(exp(-T*d))/n

% plot the final sequence and some analytics
figure()
ax = subplot(2, 1, 1);
self.plot_seq(ax);

ax = subplot(2, 1, 2);
self.plot_stats(ax);
end
