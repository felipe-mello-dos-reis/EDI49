function SecaoCriticaIntervalo

    TI = cputime;
    
    % Chama a função de entrada de dados
    [Xc, Yc, INC, Nc, gamma_c, sigma_min, sigma_max, dp, eta, F] = EntradaDeDadosSecaoCritica;

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
    coef_ang_e_min = zeros(2,length(eta));
    coef_ang_e_max = zeros(2,length(eta));
    for i = 1:length(eta)
        [M_cric(1,i), z_cric(1,i)] = max(M(:,i));
        fprintf('Fase %d:\n', i);
        fprintf('\tMomento fletor máximo: %f kN.cm\n', M_cric(1,i));
        fprintf('\tPosição do momento fletor máximo: %f cm\n', Lspan(z_cric(1,i)));

        coef_ang_e_min(1,i) = (M_cric(1,i) + Wb*sigma_min(i))/eta(i);
        coef_ang_e_min(2,i) = (M_cric(1,i) - Wt*sigma_max(i))/eta(i);
        coef_ang_e_max(1,i) = (M_cric(1,i) + Wb*sigma_max(i))/eta(i);
        coef_ang_e_max(2,i) = (M_cric(1,i) - Wt*sigma_min(i))/eta(i);
    end

    % Define the range of x and y
    invF = linspace(0,1.5*1e-3, 1001);
    e = linspace(-yt,yb,1001);
   
    

    a_e_max(1,1) = min(coef_ang_e_max(1,:));
    a_e_max(2,1) = min(coef_ang_e_max(2,:));
    a_e_max(3,1) = 0;
    a_e_min(1,1) = max(coef_ang_e_min(1,:));
    a_e_min(2,1) = max(coef_ang_e_min(2,:));
    a_e_min(3,1) = 0;
    b_e_min(1,1) = -kt;
    b_e_min(2,1) = +kb;
    b_e_min(3,1) = -(yt-dp);
    b_e_max(1,1) = -kt;
    b_e_max(2,1) = +kb;
    b_e_max(3,1) = yb-dp;
    


    fprintf('\n');
    
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


    figure(2)
    hold on
    % Create a grid of points
    [INVF, E] = meshgrid(invF, e);

    % Define the inequalities
    I = E >= a_e_min(1,1)*INVF + b_e_min(1,1);
    II = E <= a_e_max(1,1)*INVF + b_e_max(1,1);
    III = E <= a_e_max(2,1)*INVF + b_e_max(2,1);
    IV = E >= a_e_min(2,1)*INVF + b_e_min(2,1);
    Vb = E <= a_e_max(3,1)*INVF + b_e_max(3,1);
    Vt = E >= a_e_min(3,1)*INVF + b_e_min(3,1);

    % Combine the inequalities
    Z = I & II & III & IV & Vb & Vt;

    % Plot the region that satisfies all inequalities
    colormap([1 1 1; 1 1 0]);  % Set the colormap so that 0 is white and 1 is yellow
    image(invF, e, Z);
    axis xy;
    set(gca, 'Ydir', 'reverse')  % This line inverts the y-axis
    
    for i = 1:(length(eta)+1)
        plot(invF, a_e_min(i,1)*invF + b_e_min(i,1))
        plot(invF, a_e_max(i,1)*invF + b_e_max(i,1))
    end
    xlabel('1/F (kN^-1)')
    ylabel('e (cm)')
    legend('I','II','IV','III','Vt','Vb')
    title('Diagrama de Magnel')
    xlim([0 max(invF)])
    ylim([-yt yb])
    set(gca, 'Ydir', 'reverse')  % This line inverts the y-axis
    hold off;
    
    % print -depsc2 ../../images/exame_q5_a.eps
    
% Find the indices where Z is true
indices = find(Z);

% Get the corresponding values of INVF
INVF_values = INVF(indices);

% Find the maximum value
max_INVF = max(INVF_values);
% Get the corresponding values of E
E_values = E(indices);

% Find the E value corresponding to the maximum INVF value
max_E_index = find(INVF_values == max_INVF);
max_E = E_values(max_E_index);

fprintf('O valor máximo de F é %f kN/cm^2\n', 1/max_INVF);
fprintf('O valor máximo de e é %f cm\n', max_E);

    fprintf('SOLUÇÃO CONCLUÍDA\n\n');
    TF = cputime;
    fprintf('TRABALHO COMPUTACIONAL, SEGUNDOS %e\n', TF - TI);

end
