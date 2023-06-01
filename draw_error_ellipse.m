function [a, b] = draw_error_ellipse(data, sigma)
    cols = size(data,2);
    mu = mean(data,1);
    X0 = data - mu;

    conf = 2*normcdf(sigma)-1;
    scale = chi2inv(conf,cols);
    Cov = X0'*X0./size(X0,1) .* scale;

    [V, D] = eig(Cov);
    [D, order] = sort(diag(D), 'descend');
    D = diag(D);
    a = 2*sqrt(D(1,1));
    b = 2*sqrt(D(2,2));
    V = V(:,order);

    N = 100;
    t = linspace(0, 2*pi, N);
    e = [cos(t); sin(t)];
    VV = V*sqrt(D);
    e = VV*e + mu';
    
    plot(e(1,:), e(2,:), 'black');
end