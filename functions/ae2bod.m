function res = ae2bod(alpha, beta)
    % Returns the Direction Cosine Matrix that rotates from Aerodynamic
    % to Body frame.
    %
    % Inputs:
    % * alpha: Angle of attack [rad]
    % * beta:  Angle of sideslip [rad]
    %
    % Output:
    % * res:   AE2BOD Direction Cosine Matrix
    %

    sa = sin(alpha);
    ca = cos(alpha);
    sb = sin(beta);
    cb = cos(beta);

    res = [ca*cb, -ca*sb, -sa;
              sb,     cb,   0;
           sa*cb, -sa*sb,  ca];
end