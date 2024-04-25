function DisposicaoLongitudinal

    TI = cputime;
    
    % Chama a função de entrada de dados
    [Xc, Yc, INC, Nc, gamma_c, sigma_min, sigma_max, dp, eta, F] = EntradaDeDadosDisposicaoLongitudinal;

    % Traduz o sistema de coordenadas
    % Nc = size(Xc,2);

    AREA = 0;
    Sx = 0;
    Sy = 0;
    
    for I = 1:(Nc)
        ai = Xc(I) * Yc(I+1) - Xc(I+1) * Yc(I);
        AREA = AREA + ai;
        Sx = Sx + ai * (Yc(I) + Yc(I+1));
        Sy = Sy + ai * (Xc(I) + Xc(I+1));
    end

    AREA = (1 / 2) * AREA;
    Sx = (1 / 6) * Sx;
    Sy = (1 / 6) * Sy;
    
    Xcg = Sy / AREA;
    Ycg = Sx / AREA;

    % Translação das coordenadas
    for I = 1:(Nc+1)
        Xc(I) = Xc(I) - Xcg;
        Yc(I) = Yc(I) - Ycg;
    end

    AREA = 0;
    Sx = 0;
    Sy = 0;
    Ixx = 0;
    Iyy = 0;
    Ixy = 0;
    
    for I = 1:Nc
        ai = Xc(I) * Yc(I+1) - Xc(I+1) * Yc(I);
        AREA = AREA + ai;
        Sx = Sx + ai * (Yc(I) + Yc(I+1));
        Sy = Sy + ai * (Xc(I) + Xc(I+1));
        Ixx = Ixx + ai * (Yc(I)^2 + Yc(I) * Yc(I+1) + Yc(I+1)^2);
        Iyy = Iyy + ai * (Xc(I)^2 + Xc(I) * Xc(I+1) + Xc(I+1)^2);
        Ixy = Ixy + ai * (Xc(I) * Yc(I+1) + 2 * (Xc(I) * Yc(I) + Xc(I+1) * Yc(I+1)) + Xc(I+1) * Yc(I));
    end

    SINAL_DA_CIRCUICAO = AREA;

    if SINAL_DA_CIRCUICAO > 0
        AREA = (1 / 2) * AREA;
        Sx = (1 / 6) * Sx;
        Sy = (1 / 6) * Sy;
        Ixx = (1 / 12) * Ixx;
        Iyy = (1 / 12) * Iyy;
        Ixy = (1 / 24) * Ixy;
    else
        AREA = -(1 / 2) * AREA;
        Sx = -(1 / 6) * Sx;
        Sy = -(1 / 6) * Sy;
        Ixx = -(1 / 12) * Ixx;
        Iyy = -(1 / 12) * Iyy;
        Ixy = -(1 / 24) * Ixy;
    end

    yb = abs(min(Yc));
    yt = abs(max(Yc));
    I = Ixx;
    Wb = I/yb;
    Wt = I/yt;
    kb = Wt/AREA;
    kt = Wb/AREA;
    % Tr_I = Ixx + Iyy;

    % Xmax = Xc(1, 1);
    % Xmin = Xc(1, 1);

    % for I = 1:Nc
    %     if Xc(I) > Xmax
    %         Xmax = Xc(I);
    %     end

    %     if Xc(I) < Xmin
    %         Xmin = Yc(I);
    %     end
    % end

    % b = Xmax - Xmin;

    % Ymax = Yc(1, 1);
    % % Ymin = Yc(1, 1);

    % for I = 1:Nc
    %     if Yc(I) > Ymax
    %         Ymax = Yc(I);
    %     end

    %     if Yc(I) < Ymin
    %         Ymin = Yc(I);
    %     end
    % end

    % h = Ymax - Ymin;
    % Ysmin = Ymin;
    % Ysmax = Ymax;

    % if Ns > 0
    %     Ysmin = min(Ys);
    %     Ysmax = max(Ys);
    % end

    % Realiza a análise selecionada
    
    L = 2200; % cm
    q_load = 0.065; % kN/cm2
    Lspan = linspace(0,L,L+1);
    g = zeros(length(Lspan),length(eta));
    q = zeros(length(Lspan),length(eta));
    M = zeros(length(Lspan),length(eta));
    for z = 1:length(Lspan)
        g(z,1) = gamma_c*AREA;
        g(z,2) = gamma_c*AREA;
        q(z,1) = 0;
        q(z,2) = q_load;
        M(z,1) = abs((g(z,1)+q(z,1))*Lspan(z)^2/2 - gamma_c*AREA*L/2*Lspan(z));
        M(z,2) = abs((g(z,2)+q(z,2))*Lspan(z)^2/2 - (gamma_c*AREA+q_load)*L/2*Lspan(z));
    end

    M_cric = zeros(1,length(eta));
    z_cric = zeros(1,length(eta));
    e_min = zeros(length(Lspan),length(eta)*2+1);
    e_max = zeros(length(Lspan),length(eta)*2+1);

    e_min_cric = zeros(1,length(eta)*2+1);
    e_max_cric = zeros(1,length(eta)*2+1);

    for i = 1:length(eta)
        [M_cric(1,i), z_cric(1,i)] = max(M(:,i));
        fprintf('Fase %d:\n', i);
        fprintf('\tMomento fletor máximo: %f kN.cm\n', M_cric(1,i));
        fprintf('\tPosição do momento fletor máximo: %f cm\n', Lspan(z_cric(1,i)));

        e_min(:,(i-1)*2 + 1) = 1/F(1,i)*(M(:,i) + Wb*sigma_min(i)) - kt;
        e_min(:,(i-1)*2 + 2) = 1/F(1,i)*(M(:,i) - Wt*sigma_max(i)) + kb;
        e_max(:,(i-1)*2 + 1) = 1/F(1,i)*(M(:,i) + Wb*sigma_max(i)) - kt;
        e_max(:,(i-1)*2 + 2) = 1/F(1,i)*(M(:,i) - Wt*sigma_min(i)) + kb;
        e_min_cric(1,(i-1)*2 + 1) = 1/F(1,i)*(M_cric(1,i) + Wb*sigma_min(i)) - kt;
        e_min_cric(1,(i-1)*2 + 2) = 1/F(1,i)*(M_cric(1,i) - Wt*sigma_max(i)) + kb;
        e_max_cric(1,(i-1)*2 + 1) = 1/F(1,i)*(M_cric(1,i) + Wb*sigma_max(i)) - kt;
        e_max_cric(1,(i-1)*2 + 2) = 1/F(1,i)*(M_cric(1,i) - Wt*sigma_min(i)) + kb;
    end

    e_min(:,length(eta)*2+1) = -(yt-dp);
    e_max(:,length(eta)*2+1) = yb-dp;

    e_min_z = max(e_min, [], 2);  % Find the maximum of each row in e_min
    e_max_z = min(e_max, [], 2);  % Find the minimum of each row in e_max
    

    
    
    e_min_cric(1,length(eta)*2+1) = -(yt-dp);
    e_max_cric(1,length(eta)*2+1) = yb-dp;
    
    e_min_f = max(e_min_cric);
    e_max_f = min(e_max_cric);
    fprintf('\n');


    fprintf('e >= %f cm\n', e_min_f);
    fprintf('e <= %f cm\n', e_max_f);

    for I = 1:5
        fprintf('\n');
    end
    
    close all
    figure(1)
    plot(Xc, Yc, '-b')
    xlabel('x (cm)')
    ylabel('y (cm)')
    title('Secao Transversal')
    axis('equal')
    xlim([min(min(Xc, Yc)) max(max(Xc, Yc))]*1.25)
    ylim([min(min(Xc, Yc)) max(max(Xc, Yc))]*1.25)
    % print -depsc2 ../../images/exame_q5_a.eps


    figure(2)
    plot(Lspan, -e_min_z, '-r', Lspan, -e_max_z, '-b', Lspan, yt*ones(size(Lspan)), '-k', Lspan, -yb*ones(size(Lspan)), '-k', Lspan, zeros(size(Lspan)), '--k')
    xlabel('x (cm)')
    ylabel('e (cm)')
    title('Região Limite')
    legend('e_{min}', 'e_{max}','Borda superior: y=yt', 'Borda inferior: y=-yb', 'CG')
    xlim([0 L])
    ylim([-yb*1.25 yt*1.25])
    grid on


    fprintf('SOLUÇÃO CONCLUÍDA\n\n');
    TF = cputime;
    fprintf('TRABALHO COMPUTACIONAL, SEGUNDOS %e\n', TF - TI);

end
