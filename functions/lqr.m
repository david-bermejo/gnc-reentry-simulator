function [K, P] = lqr(A, B, Q, R)
    % Computes the feedback matrix for the continuous time, infinite-
    % horizon Linear Quadratic Regulator.
    %
    % The gain matrix is: K = inv(R)*B'*P
    %
    % where P is the solution of the following algebraic Ricatti equation:
    % 0 = P*A + A'*P -P*B*inv(R)*B'*P + Q
    %
    % NOTE: this function is compatible with MATLAB's code generator.
    %

    rows = size(B,2);
    cols = size(A,1);
    K = zeros(rows,cols);

    J = R\B';
    G = B*J;

    H = [A -G; -Q -A'];
    [V, D] = eig(H);
    idx = find(diag(D) <= 0);
    sz = size(D,1);
    
    U11 = V(1:(sz/2), idx);
    U21 = V((sz/2+1):sz, idx);
    
    P = real(U21/U11);
    K(:,:) = J*P;
end