syms x y z xk yk zk vx vy vz vxk vyk vzk mu qw qx qy qz

qwx = qw*qx; qwy = qw*qy; qwz = qw*qz;
qxx = qx^2; qxy = qx*qy; qxz = qx*qz;
qyy = qy^2; qyz = qy*qz;
qzz = qz^2;

dcm = [1-2*(qyy+qzz)    2*(qxy - qwz)    2*(qwy + qxz);
       2*(qxy + qwz)    1-2*(qxx+qzz)    2*(qyz - qwx);
       2*(qxz - qwy)    2*(qwx + qyz)    1-2*(qxx+qyy)];

g = -mu/(x^2+y^2+z^2)^(3/2) * [x; y; z];
vp = dcm * [x; y; z];

Jg = simplify(jacobian(g, [x y z]));
Jqrot = simplify(jacobian(vp, [qw qx qy qz]));

r = [x; y; z];
rk = [xk; yk; zk];
v = [vx; vy; vz];

Rk = sqrt((x-xk)^2 + (y-yk)^2 + (z-zk)^2);
%Vk = (r - rk)'*v/Rk;
Vk = ((x-xk)*(vx-vxk) + (y-yk)*(vy-vyk) + (z-zk)*(vz-vzk))/sqrt((x-xk)^2 + (y-yk)^2 + (z-zk)^2);

JRk = simplify(jacobian(Rk, [x y z]));
JVk = simplify(jacobian(Vk, [x y z vx vy vz]));
JRk
JVk