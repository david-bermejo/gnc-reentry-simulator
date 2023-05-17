function Cb = ctrl(A, B)
    % Computes the controlability matrix of a given system in its
    % state space model form:
    %
    %   x_dot = A*x + B*u
    %       y = C*x + D*u
    %
    % The controllability matrix is calculated as follows:
    % Cb = [B A*B A^2*B ... A^(n-1)*B]
    %
    % NOTE: This function is compatible with MATLAB's code generator.

    sz = size(A,1);
    Cb = zeros(sz,sz);
    Cb(:,1) = B;

    for i=2:sz
        Cb(:,i) = A*Cb(:,i-1);
    end
end