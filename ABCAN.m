%%% This is an implementation of ABC algorithm in order to search some
%%% parameters in a dynamic system to simulate a signal.

%%% Inputs
%%% SN = Number of population (single value)
%%% limit = Number of limit (single value)
%%% MFE = Number of generations (single value)
%%% N_exp = Number of iterations (single value)
%%% coupled = is attached with windmi system (and with the number of
%%% parameters to search in the system)  (single value)
%%% n_best = Number of solutions selected (single value)

%%% Outputs
%%% XI = All signal solutions given in a matrix (N_exp,MFE,size of the solution vector)
%%% PI = All parameters solutions given in a matrix (N_exp,MFE,size of the parameters vector)
%%% EI = Value of the objective function in the last generation
%%% best_Xs = Best signal solutions given in a matrix (n_best,size of the solution vector)
%%% best_Ps = Best parameters solutions given in a matrix (n_best,size of the parameters vector)


function [XI,PI,EI,best_Xs,best_Ps]=ABCAN(SN,limit,MFE,N_exp,coupled,n_best)
    rng(7)
    if coupled
        n_params = 16;
    else    
        n_params = 7;
    end
%     bounds = [0 2; 0.1 0.5; 10 20; 15 23; 0 1;  0 1; 0 1; 0 2; 0.1 0.5; 10 20; 15 23];
    bounds = [0 2; 0.1 0.5; 10 20; 15 23; 0 1;  0 1; 0 1];
    n_ex = 0;
    EI = zeros(N_exp, MFE);
    PI = zeros(N_exp, MFE, n_params);
    XI = zeros(N_exp, MFE, 5478);
    tic;
    while(n_ex~=N_exp)%%% Número de iteraciones
        g = 0;
        P = Inicializacion(SN, n_params, bounds);
        X = windmi(P, coupled);
        FIT = FunObjA(X);
        tic;
        while(g~=MFE)
            for s=1:SN
              u = (rand*2)-1;
              idxs = 1:SN;
              idxs = idxs(idxs ~= s);
              i = 1 + round(rand*(SN-2));
              i = idxs(i);
              mut = mutacion(P(s,:), P(i,:), u);
              V = windmi(mut, coupled);
              EV = FunObjA(V);
              if EV < FIT(s)
                P(s,:) = mut;
                X(s,:) = V;
                FIT(s) = EV;
              end
            end
            cont = zeros(1,SN);
            Prob = FIT ./ sum(FIT);

            for s=1:SN
                prob_cross = rand();
                if Prob(s) > prob_cross
                      u = (rand*2)-1;
                      idxs = 1:SN;
                      idxs = idxs(idxs ~= s);
                      i = 1 + round(rand*(SN-2));
                      i = idxs(i);
                      son = mutacion(P(s,:), P(i,:), u);
                      V1 = windmi(son, coupled);
                      EV1 = FunObjA(V1);
                      if EV1 < FIT(s)
                        P(s, :) = son;
                        FIT(s) = EV1;
                        X(s, :) = V1;
                      end
                else
                     cont(s)=cont(s)+1;
                     if cont(s) >= limit
                         P(s, :) = Inicializacion(1, n_params, bounds);
                         X(s, :) = windmi(P(s, :), coupled);
                         FIT(s) = FunObjA(X(s, :));
                     end
                end
            end
            [m,I]=min(FIT);
            EI(n_ex + 1, g + 1) = m;
            XI(n_ex + 1, g + 1, :) = X(I, :);
            PI(n_ex + 1, g + 1, :) = P(I,:);
            g = g + 1;
            disp(['Experimento: ' num2str(n_ex + 1) ', Generación: ' num2str(g) ', MinError: ' num2str(m)])
        end
        n_ex = n_ex + 1;
    end
    figure(1)
    plot(EI')
    title(['Función de error para los ' num2str(N_exp) ' experimentos'])
    best_Xs = zeros(n_best, 5478);
    best_Ps = zeros(n_best, n_params);
    EIt = EI;
    for i=1:n_best
       [best_exps, best_exp_idxs] = min(EIt,[],2);
       [~, best_exp] = min(best_exps);
       best_exp_idx = [best_exp best_exp_idxs(best_exp)];
       best_params = squeeze(PI(best_exp_idx(1),best_exp_idx(2),:));
       best_series = squeeze(XI(best_exp_idx(1),best_exp_idx(2),:));
       EIt(EIt == EIt(best_exp_idx(1),best_exp_idx(2))) = nan;
       best_Xs(i,:) = best_series;
       best_Ps(i,:) = best_params;
    end
    figure(2)
    plot(best_Xs(1,:))
    hold on
    load('señal.mat','ft');
    plot(ft/max(ft))
    hold off
    
    save('XI.mat','XI');
    save('PI.mat','PI');
    save('EI.mat','EI');
    save('best_series.mat','best_Xs');
    save('best_params.mat','best_Ps');
    toc;
end