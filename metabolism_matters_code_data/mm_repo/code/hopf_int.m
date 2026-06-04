function [FC,CV,Cvth,A]=hopf_int(gC,f_diff,sigma,a_offset)
% HOPF_INT  Linearised Stuart-Landau / Hopf solution for FC / COV / Jacobian.
%
%   sigma : scalar (homogeneous noise) OR length-N vector (per-region noise).
%   a_offset (optional) : length-N vector added element-wise to the scalar
%                         bifurcation parameter a = -0.02, i.e. a_i = a + a_offset(i).
%                         Used by the NEMO / membrane-polarization perturbation
%                         in Ceff_Computation.m. Defaults to zeros.

    % --- DEBUG CODE ---
    if any(f_diff == 0) || isempty(f_diff)
        disp('CRITICAL: f_diff passed to hopf_int contains zeros or is empty!');
        disp(['Size of f_diff: ', num2str(size(f_diff))]);
    end
    % ------------------

    N=size(gC,1);
    a=-0.02;

    if nargin < 4 || isempty(a_offset)
        a_offset = zeros(N,1);
    end
    a_vec = a + a_offset(:);                   % N x 1 per-region bifurcation

    wo = f_diff'*(2*pi);

    Cvth = zeros(2*N);

    % Jacobian:

    s = sum(gC,2);
    B = diag(s);

    Axx = diag(a_vec) - B + gC;                % was a*eye(N)
    Ayy = Axx;
    Axy = -diag(wo);
    Ayx = diag(wo);

    A = [Axx Axy; Ayx Ayy];

    % Noise covariance: scalar sigma -> uniform; vector sigma -> per-region
    if isscalar(sigma)
        Qn = (sigma^2) * eye(2*N);
    else
        sv = sigma(:);
        if numel(sv) ~= N
            error('hopf_int: sigma must be scalar or length-N vector (got %d for N=%d).', numel(sv), N);
        end
        Qn = diag([sv.^2; sv.^2]);
    end
    
    % Check stability of the origin:
    [~,d] = eig(A);
    d = diag(d);
    Remax = max(real(d));
    if Remax >= 0
      disp('Warning: the origin is not stable') 
      %FC = []; CV=[]; Cvth=[]; A=[];
      return
end

%Cvth=sylvester(A,A',-Qn);
Cvth = lyap(A, Qn);
FCth=corrcov(Cvth);
FC=FCth(1:N,1:N);
CV=Cvth(1:N,1:N);
