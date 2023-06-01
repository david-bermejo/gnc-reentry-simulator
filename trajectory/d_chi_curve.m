function d_chi = d_chi_curve(d_chi0, d_chif, M0, Mf, M)
    if M > M0
        d_chi = d_chi0;
    elseif M < Mf
        d_chi = d_chif;
    else
        d_chi = d_chi0 + (d_chif - d_chi0)/(Mf - M0)*(M - M0);
    end
end