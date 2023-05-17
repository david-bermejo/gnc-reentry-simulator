function dcm = quat2dcm(q)
    % Transforms a quaternion into its equivalent Direction Cosine Matrix.
    
    qwx = q(1)*q(2); qwy = q(1)*q(3); qwz = q(1)*q(4);
    qxx = q(2)^2;    qxy = q(2)*q(3); qxz = q(2)*q(4);
    qyy = q(3)^2;    qyz = q(3)*q(4);
    qzz = q(4)^2;
    
    dcm = [1-2*(qyy+qzz)    2*(qxy - qwz)    2*(qwy + qxz);
           2*(qxy + qwz)    1-2*(qxx+qzz)    2*(qyz - qwx);
           2*(qxz - qwy)    2*(qwx + qyz)    1-2*(qxx+qyy)];
end