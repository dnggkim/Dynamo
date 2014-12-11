% Example 0 for DPG 2012 Stuttgart

randseed(2783);

% 2-qubit Heisenberg chain, XY controls at one end, state transfer |00> -> |11>

q = 2;
dim = 2 * ones(1, q);

desc = sprintf('Isotropic Heisenberg chain, %d qubits, XY control at one end.', q);
fprintf('%s\n\n', desc);

J = 2 * [1 1 1];
C = diag(ones(1, q-1), 1); % topology: linear chain
H = heisenberg(dim, @(s,a,b) J(s)*C(a,b));
[C, cl] = control_ops(dim, '1xy');

initial = [1 0 0 0].';
final = [0 0 0 1].';

dyn = dynamo('closed ket', initial, final, H, C);
dyn.system.set_labels(desc, dim, cl);
dyn.seq_init(100, 6 * [1, 0]);
dyn.easy_control([-0.1, 0.05]);

dyn.ui_open();
dyn.search_BFGS();
%dyn.analyze();
%figure; dyn.plot_X();
%figure; dyn.plot_seq();
