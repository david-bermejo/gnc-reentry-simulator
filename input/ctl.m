function out = ctl()
    out.freq = 20;
    out.tsamp = 1/out.freq;

    %% Longitudinal Controller Weight Matrices
    out.Q_lon = zeros(2);
    out.R_lon = zeros(2);

    out.Q_lon(1,1) = 1/deg2rad(5)^2;
    out.Q_lon(2,2) = 1/deg2rad(0.1)^2;

    out.R_lon(1,1) = 1/40^2;
    out.R_lon(2,2) = 1/10400^2;

    %% Lateral Controller Weight Matrices
    out.Q_lat = zeros(4);
    out.R_lat = zeros(4);

    out.Q_lat(1,1) = 1/deg2rad(5)^2;
    out.Q_lat(2,2) = 1/deg2rad(5)^2;
    out.Q_lat(3,3) = 1/deg2rad(0.1)^2;
    out.Q_lat(4,4) = 1/deg2rad(0.1)^2;

    out.R_lat(1,1) = 1/40^2;
    out.R_lat(2,2) = 1/40^2;
    out.R_lat(3,3) = 1/1600^2;
    out.R_lat(4,4) = 1/7600^2;
end