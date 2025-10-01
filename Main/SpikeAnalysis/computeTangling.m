function Q = computeTangling(X, epsilon)
% X: latent_dim x time x trials
% epsilon: small constant (e.g., 1e-6)
% Q: 1 x trials vector, Q(i) is the maximal tangle for trial i
% We quantified trajectory tangling using QðtÞ = maxt0 kx_t  x_t0 k 2 kxt
% xt0 k 2 + ε ; (Equation 1) where xt is the neural state at time t (i.e.,
% a vector containing the neural responses at that time), x_t is the
% temporal derivative of the neural state, k, k is the Euclidean norm, and
% ε is a small constant that prevents division by zero (STAR Methods). QðtÞ
% becomes high if there exists a state at a different time, t 0 , that is
% similar but associated with a dissimilar derivative. We take the maximum
% to ask whether the state at time t ever becomes tangled with any other
% state. This maximum is taken with t 0 indexing across time during all
% conditions.
if nargin < 2
    epsilon = 1e-6;
end

[latent_dim, T, nTrials] = size(X);
Q = zeros(T, nTrials);

for trial = 1:nTrials
    X_trial = X(:, :, trial);  % [latent_dim x time]
    dX = diff(X_trial, 1, 2);  % [latent_dim x T-1]
    dX = [dX, dX(:,end)];      % Pad to length T

    Qt = zeros(1, T);

    for t = 1:T
        xt = X_trial(:, t);
        dxt = dX(:, t);

        xt_xtp = xt-X_trial;        % [latent_dim x T]
        dxt_dxtp = dxt-dX;          % [latent_dim x T]

        num = sum(dxt_dxtp .^ 2, 1);     % Numerator: diff of derivatives
        denom = sum(xt_xtp .^ 2, 1) + epsilon;  % Denominator: state dist + eps

%         num(t) = -inf;  % Don't compare to self

        Qt(t) = max(num ./ denom);    % For timepoint t, maximal tangle
    end

    Q(:,trial) = Qt;   % For this trial, report maximal Q across all t
end

end
