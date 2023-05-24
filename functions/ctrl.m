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

    sz1 = size(A,1);
    sz2 = size(B,2);
    Cb = zeros(sz1,sz1*sz2);
    Cb(:,1:sz2) = B;

    for i=2:sz1
        idx = (i-1)*sz2;
        pidx = (i-2)*sz2;
        Cb(:,idx+1:idx+sz2) = A*Cb(:,pidx+1:pidx+sz2);
    end
end