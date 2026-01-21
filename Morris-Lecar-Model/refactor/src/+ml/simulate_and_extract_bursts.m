function raw = simulate_and_extract_bursts(I_values, baseParams, sim)
%SIMULATE_AND_EXTRACT_BURSTS Simulate Morris–Lecar and extract one burst per I.
%
% raw = struct with fields:
%   .I_values  (Nx1)
%   .periods   (Nx1)
%   .spikeCounts (Nx1)
%   .segT      (cell Nx1) time of extracted burst segment
%   .segZ      (cell Nx1) [V w y] segment states
%   .fullT     (cell Nx1) full simulation time (optional)
%   .fullZ     (cell Nx1) full simulation states (optional)
%

N = numel(I_values);
raw = struct();
raw.I_values = I_values(:);
raw.periods = zeros(N,1);
raw.spikeCounts = zeros(N,1);
raw.segT = cell(N,1);
raw.segZ = cell(N,1);
raw.fullT = cell(N,1);
raw.fullZ = cell(N,1);

% ODE options: keep defaults; you can expose this via sim.odeopts later.
%odeopts = odeset('RelTol',1e-6,'AbsTol',1e-8);

for k = 1:N
    p = baseParams; p.I = I_values(k);

    % Solve ODE RHS
    [t, Z] = ode15s(@(t,z) ml.rhs(t, z, p), sim.tspan, sim.z0);
    raw.fullT{k} = t;
    raw.fullZ{k} = Z;

    % Remove transient
    keepIdx = t > sim.Tfinal*sim.transientFraction;
    t2 = t(keepIdx);
    Z2 = Z(keepIdx,:);

    % Extract last burst period usind threshold detection
    [segT, segZ] = ml.extract_last_burst(t2, Z2);
    raw.segT{k} = segT;
    raw.segZ{k} = segZ;
    raw.periods(k) = segT(end) - segT(1);

    % Test spike count consistency
    raw.spikeCounts(k) = ml.count_spikes(segZ(:,1));

    fprintf('I=%.4f : burst dur=%.3f, spikes=%d\n', p.I, raw.periods(k), raw.spikeCounts(k));
end

end
