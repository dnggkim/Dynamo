% Example 2 for DPG 2012 Stuttgart

randseed(2783);

% 3-qubit Heisenberg chain, XY controls at one end, QFT gate
q = 3;
dim = 2 * ones(1, q);

desc = sprintf('Isotropic Heisenberg chain, %d qubits, XY control at one end.', q);
fprintf('%s\n\n', desc);

J = 2 * [1 1 1]; % Heisenberg interaction
C = diag(ones(1, q-1), 1); % topology: linear chain
H = heisenberg(dim, @(s,a,b) J(s)*C(a,b));
[C, cl] = control_ops(dim, '1xy');

final = qft(q);
initial = eye(size(final));

dyn = dynamo('S gate', initial, final, H, C);
dyn.system.set_labels(desc, dim, cl);
dyn.seq_init(100, 16 * [1, 0]);
dyn.easy_control(0.1 * ones(1,2));

dyn.ui_open();
dyn.search_BFGS(dyn.full_mask(), struct('Display', 'final', 'plot_interval', 1));
%dyn.analyze();
