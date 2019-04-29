%%% This is the implementation of two different objective functions

%%% Inputs
%%% X = Vector with the signal given by the dynamic system

%%% Outputs
%%% FIT = Retutn a single value of the objectie function


function [FIT] = FunObjA(X)

    load ('señal.mat','ft');
    ft = ft/max(ft);
    
    MSE = 0.5 * mean((ft' - X).^2,2);
    NSE = 1 - (sum((ft' - X).^2,2)./(sum((ft'- mean(ft')).^2)));
    MAE = (1/5478)*(sum((abs(ft' - X)),2));
    diffft = diff(ft');
    diffX = diff(X,1,2);
    BK = (sign(diffft) == sign(diffX));
    POCID = 1/5478 * sum(BK, 2);
    neg_pen = sum(X < 0, 2);
    
    %%% Objective functions
    % 1
    FIT = ((1-NSE).*MAE) ./ (POCID);
    %2
%     FIT = (((1-NSE).*MAE+ neg_pen + MSE) ./ (POCID)); 
    %%%%%%
end
%%%%%%%%%%%
    
