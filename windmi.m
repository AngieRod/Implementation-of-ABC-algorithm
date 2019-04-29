%%% This is an implementation of the Multiscroll caotic system based in
%%% Chua chaotic circuit - also contain an coupled system

%%% Inputs
%%% P = vector of parameters
%%% Coupled = Single value (1 or 0), indicating if is coupled or not.

%%% Outputs
%%% X, Y, Z = State variables of the system. X represent the main signal.

function [X,Y,Z]=windmi(P, coupled)
    N = 5478;
%     N = 6100;
    %%%
    XA = zeros(size(P,1),N);
    YA = zeros(size(P,1),N);
    ZA = zeros(size(P,1),N);
    
    %%%
    a = P(:,1);
    b = P(:,2);
    alpha = P(:,3);
    beta = P(:,4);
    x0 = P(:,5);
    y0 = P(:,6);
    z0 = P(:,7);
    
    XA(:,1) = x0;
    YA(:,1) = y0;
    ZA(:,1) = z0;
    
    if coupled 
        XU = zeros(size(P,1),N);
        YU = zeros(size(P,1),N);
        ZU = zeros(size(P,1),N);
        a1 = P(:,8);
        b1 = P(:,9);
        alpha1 = P(:,10);
        beta1 = P(:,11);
        x01 = P(:,12);
        y01 = P(:,13);
        z01 = P(:,14);
        kx = P(:,15);
        k1x = P(:,16);
        
        XU(:,1) = x01;
        YU(:,1) = y01;
        ZU(:,1) = z01;                
    end
    %%%    
  
    T = 0.1;
    for i=1:N-1
        if coupled            
            XA(:,i+1) = min(max(alpha.*YA(:,i)*T + alpha.*b.*T.*sin(pi*XA(:,i)./(2*a)) + XA(:,i),-1e10),1e10);
            YA(:,i+1) = min(max(XA(:,i).*T + ZA(:,i).*T - YA(:,i).*T + YA(:,i)+ kx .* (YU(:,i)),-1e10),1e10);
            ZA(:,i+1) = min(max(-beta.*T.*YA(:,i) + ZA(:,i),-1e10),1e10);
            
            XU(:,i+1) = min(max(alpha1.*YU(:,i)*T + alpha1.*b1.*T.*sin(pi*XU(:,i)./(2*a1)) + XU(:,i),-1e10),1e10);
            YU(:,i+1) = min(max(XU(:,i).*T + ZU(:,i).*T - YU(:,i).*T + YU(:,i)+ k1x .* (YA(:,i)),-1e10),1e10);
            ZU(:,i+1) = min(max(-beta1.*T.*YU(:,i) + ZU(:,i),-1e10),1e10);
            
        else
            XA(:,i+1) = min(max(alpha.*YA(:,i)*T + alpha.*b.*T.*sin(pi*XA(:,i)./(2*a)) + XA(:,i),-1e6),1e6);
            YA(:,i+1) = min(max(XA(:,i).*T + ZA(:,i).*T - YA(:,i).*T + YA(:,i),-1e6),1e6);
            ZA(:,i+1) = min(max(-beta.*T.*YA(:,i) + ZA(:,i),-1e6),1e6);
        end
    end
    X = XA;
    Y = YA;
    Z = ZA;
    
    X = X ./max(abs(X(1:5478)), [], 2);
    Y = Y ./max(abs(Y(1:5478)), [], 2);
    Z = Z ./max(abs(Z(1:5478)), [], 2);
    
end