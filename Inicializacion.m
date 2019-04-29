function [Pob]=Inicializacion(SN, n_params, bnds)
    bounds = zeros(n_params, 2);
    if size(bnds, 1) < n_params
        bounds(1:size(bnds, 1), :) = bnds;
        bounds(size(bnds, 1) + 1:end, :) = repmat([0 1], [n_params-size(bnds, 1), 1]);
    elseif size(bnds, 1) == n_params
        bounds = bnds;
    else
        disp('NOOOOOO')
    end
    rng('shuffle');
    Pob = rand(SN, n_params) .*((bounds(:, 2) - bounds(:, 1))' + bounds(:, 1)');
end