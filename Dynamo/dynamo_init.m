function dynamo_init(task, initial, final, H_drift, H_ctrl, L_drift)
% Initializes Dynamo for a system and an optimization task.
%
% Governing equation: \dot(X)(t) = -(A +\sum_k u_k(t) B_k) X(t) = -G(t) X(t)

% Ville Bergholm 2011


% Dynamo version
version = '1.3 alpha10';

% All definitions are in a global variable called OC
global OC; % and now we can access it too

if numel(OC)==0
    disp (' ');
    fprintf('DYNAMO - Quantum Dynamic Optimization Package v%s\n', version);
    disp (' ');
    disp (' ');
    disp ('(c) Shai Machnes et al. 2010-2011');
    disp ('email: shai.machnes at uni-ulm.de');
    disp (' ');
    disp ('All computer programs/code/scripts are released under the terms of the GNU Lesser General Public License 3.0 and Creative-Commons Attribution Share-Alike (see "LICENSE.txt" for details).');
    disp ('  ');
    disp ('If you use DYNAMO in your research, please add an attribution in the form of the following reference: S. Machnes et al, arXiv 1011.4874');
    disp (' ');
    disp ('For the latest version of this software, guides and information, visit http://www.qlib.info');
    disp ('  ');
    disp ('DYNAMO initialized successfully.');
    disp ('  ');    
    disp ('  ');
    drawnow;
end

task = lower(task);

%% Some basic data provenance

OC.config.version = version;
% Local time. TODO UTC or local time with timezone specifier would be better, but apparently MATLAB doesn't do that.
OC.config.date = datestr(now(), 31);
OC.config.task = task;

OC.config.expmFunc = @expm;

% TODO temporary fix: sparse to full
H_drift = full(H_drift);
for k = 1:length(H_ctrl)
    H_ctrl{k} = full(H_ctrl{k});
end

input_dim = [size(initial, 2), size(final, 2)]; % check the validity of the inputs

[system, rem] = strtok(task);
[task, rem] = strtok(rem);
[phase, rem] = strtok(rem);
out = '\nOptimize a control sequence to obtain the given';
switch system
  case {'s'}
    %% Closed system S
    if nargin ~= 5
        error('Too many parameters.')
    end
    
    switch task
      case 'state'
        out = strcat(out, ' mixed state transfer');
        % TODO more efficient Hilbert space implementation?
        OC.system = system_vec(OC.system, initial, final);
        OC.system = system_liouville(OC.system, H_drift, 0, H_ctrl);
        % g is always real, positive in this case so error_abs would work just as well
        OC.config.error_func = @error_real;
        
      case {'ket', 'gate'}
        if strcmp(task, 'ket')
            out = strcat(out, ' pure state transfer');
            if any(input_dim ~= 1)
                error('Initial and final states should be normalized kets.')
            end
        else
            out = strcat(out, ' unitary gate');
            if any(input_dim == 1)
                error('Initial and final states should be unitary operators.')
            end
        end

        OC.system.X_initial = initial;
        OC.system.X_final   = final;
        OC.system = system_hilbert(OC.system, H_drift, H_ctrl);
        
        if strcmp(phase, 'phase')
            out = strcat(out, ' (with global phase (NOTE: unphysical!))');
            OC.config.error_func = @error_real;
        else
            out = strcat(out, ' (ignoring global phase)');
            OC.config.error_func = @error_abs;
        end
        
      otherwise
        error('Unknown task.')
    end

    % global maximum of the quality function (fidelity or Q in the docs FIXME)
    OC.system.max_Q = sqrt(norm2(OC.system.X_initial) / norm2(OC.system.X_final));

    out = strcat(out, ' in a closed system.\n');

    % L: X_final' propagated backward 
    OC.cache.L_end = OC.system.X_final';

    % the generator is always Hermitian and thus normal => use exact gradient
    OC.config.gradientFunc = @gradient_exact;
    OC.config.calcPfromHfunc = @calcPfromH_exact_gradient; % When computing exact gradient, we get exponentiation for free due to the eigendecomposition (see paper for details)    

    
  case {'sb'}
    %% Open system S with bath B
    switch task
      case 'state'
        out = strcat(out, ' quantum state transfer');
        OC.system = system_vec(OC.system, initial, final);

      case 'gate'
        out = strcat(out, ' quantum map');
        if any(input_dim == 1)
            error('Initial and final states should be operators.')
        end
        
        OC.system.X_initial = lrmul(initial, initial'); % == kron(conj(initial), initial);
        OC.system.X_final   = lrmul(final, final'); % == kron(conj(final), final);
        
      otherwise
        error('Unknown task.')
    end
    out = strcat(out, ' in an open system under Markovian noise.\n');
    
    OC.system = system_liouville(OC.system, H_drift, L_drift, H_ctrl);

    % The generator isn't usually normal, so we cannot use the exact gradient method
    OC.config.error_func = @error_open;
    OC.opt.max_violation = 0;

    % L: reverse propagator
    OC.cache.L_end = eye(length(OC.system.X_final));
    
    %OC.config.error_func = @error_real; % TEST, requires also OC.cache.L{end} = X_final'
    %OC.config.gradientFunc = @gradient_first_order_aprox;
    OC.config.calcPfromHfunc = @calcPfromH_expm;

  case {'se'}
    %% Closed system S + environment E
    error('Not implemented yet.')
  
  case {'seb'}
    %% Open system S + environment E with bath B
    error('Not implemented yet.')
    
  otherwise
    error('Unknown system specification.')
end

fprintf(out);
fprintf('Optimization system dimension: %d\n', length(OC.system.X_final));

% Calculate the squared norm |X_final|^2 to scale subsequent fidelities.
% We use the Hilbert-Schmidt inner product (and the induced Frobenius norm) throughout the code.
OC.system.norm2 = norm2(OC.system.X_final);
end




function sys = system_vec(sys, initial, final)
% Set up the vec representation for the initial and final states of
% a system in Liouville space.

  % state vectors are converted to state operators
  if size(initial, 2) == 1
    initial = initial * initial';
  end
  if size(final, 2) == 1
    final = final * final';
  end

  sys.X_initial = vec(initial);
  sys.X_final   = vec(final);
end



function sys = system_hilbert(sys, H_drift, H_ctrl)
% Set up Hilbert space generators for a system.

  % (NOTE: generators are not pure Hamiltonians, there's an extra 1i!)
  sys.A = 1i * H_drift;

  n_controls = length(H_ctrl);
  sys.B = cell(1, n_controls);
  sys.B_is_superop = false(1, n_controls);
  for k=1:n_controls
      sys.B{k} = 1i * H_ctrl{k};
  end
  sys.M = inprod_B(sys.B);
end


function M = inprod_B(B)
% Computes the inner product matrix of the control operators.

  n_controls = length(B);
  M = zeros(n_controls);
  for j = 1:n_controls
    for k = 1:n_controls
      M(j,k) = inprod(B{j}, B{k});
    end
  end
  % FIXME what about dissipative controls? superoperators?
end


function sys = system_liouville(sys, H_drift, L_drift, H_ctrl)
% Set up Liouville space generators for a system.

  sys.A = -L_drift +1i*comm(H_drift);

  n_controls = length(H_ctrl);
  sys.B = cell(1, n_controls);
  sys.B_is_superop = false(1, n_controls);

  % Liouville space dimension
  dim = length(sys.X_final);

  for k=1:n_controls
      % check for Liouvillian controls
      if length(H_ctrl{k}) ~= dim
          sys.B{k} = 1i*comm(H_ctrl{k}); % Hamiltonian
      else
          sys.B{k} = -H_ctrl{k}; % Liouvillian
          sys.B_is_superop(k) = true;
      end
  end
end

