function Ob = obsv(A, C)
    % Computes the observability matrix of a given system in its
    % state space model form:
    %
    %   x_dot = A*x + B*u
    %       y = C*x + D*u
    %
    % The observability matrix is calculated as follows:
    % Ob = [C; C*A; C*A^2; ... C*A^(n-1)]
    %
    % NOTE: This function is compatible with MATLAB's code generator.

    sz1 = size(A,1);
    sz2 = size(C,1);
    Ob = zeros(sz1*sz2,sz1);

    Ob(1:sz2,:) = C;
    for i=2:sz1
        idx = (i-1)*sz2;
        pidx = (i-2)*sz2;
        Ob(idx+1:idx+sz2,:) = Ob(pidx+1:pidx+sz2,:)*A;
    end
end