function [Xc, Yc, INC, Nc, gamma_c, sigma_min, sigma_max, dp, eta, F] = EntradaDeDadosSecaoCritica

    
    %% COORDENADAS DOS PONTOS DA SEÇÃO POLIGONAL DE CONCRETO
    
    % PONTO     Xc(cm)     Yc(cm)
    
    Xc(1) = -10.0;
    Yc(1) = 0.0;
    Xc(2) = 10.0;
    Yc(2) = 0.0;
    Xc(3) = 10.0;
    Yc(3) = 82.0;
    Xc(4) = 61.0;
    Yc(4) = 94.0;
    Xc(5) = 61.0;
    Yc(5) = 103.0;
    Xc(6) = -61.0;
    Yc(6) = 103.0;
    Xc(7) = -61.0;
    Yc(7) = 94.0;
    Xc(8) = -10.0;
    Yc(8) = 82.0;
    Xc(9) = Xc(1);
    Yc(9) = Yc(1);
    

    %% INCIDENCIA DAS ARESTAS DA SEÇÃO POLIGONAL DE CONCRETO
    % ARESTA        PONTO1          PONTO2
    INC(1,1) = 1;
    INC(2,1) = 2;
    INC(1,2) = 2;
    INC(2,2) = 3;
    INC(1,3) = 3;
    INC(2,3) = 4;
    INC(1,4) = 4;
    INC(2,4) = 5;
    INC(1,5) = 5;
    INC(2,5) = 6;
    INC(1,6) = 6;
    INC(2,6) = 7;
    INC(1,7) = 7;
    INC(2,7) = 8;
    INC(1,8) = 8;
    INC(2,8) = 1;
    INC = INC';
    Nc = length(Xc)-1;
    
    %% PROPRIEDADES DO CONCRETO     
    % Efeito Rüsch   fck (kN/cm2)        gamac    
    gamma_c = 25/1e6; % kN/cm3
    sigma_min = [-2.2*1e-1, -5.0*1e-1]; % kN/cm2
    sigma_max = [16.67*1e-1, 15.0*1e-1]; % kN/cm2
    dp = 10; % cm
    eta = [1, 0.83];
    F = 2*1e3*eta; % kN/cm2


    