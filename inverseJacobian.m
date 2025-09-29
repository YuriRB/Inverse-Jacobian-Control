%%
%Aula 11: Controle no espeço cartesiano com jacobiano inverso
% Parâmetros do robô RRP (com dinâmica e controle PD)
% Parâmetros do robô RRP

% Rodar para levantar equações algébricas
% colocar tabelas com parametros
% e resultados

robot = importrobot('omronEcobra600.urdf');

% Inicialize as juntas do robô
robotConfig = homeConfiguration(robot);


% l1, l2, l3: comprimentos dos elos (m)
l1 = 0.325; l2 = 0.275; l3 = 0.210; 
% z0: altura inicial do efetor final (ponto da ferramenta)
z0 = 0.387;
% m1, m2, m3: massas dos elos (kg)
m1 = 2.0; m2 = 2.0; m3 = 0.4;
% I1, I2, I3: momentos de inércia dos elos em relação às juntas (kg·m²)
I1 = 0.01; I2 = 0.01; I3 = 0.001;
% g: aceleração da gravidade (m/s²)
g = 9.81;



% Controle PD corrige erro de posição (Kp) e erro de velocidade (Kd). Amortecimento ajuda a evitar oscilações e instabilidades.

% Ganhos do controlador PD (Proporcional-Derivativo) no espaço operacional
% Kp: matriz diagonal com ganhos proporcionais para cada junta
Kp = diag([120, 120, 120]);
% Kd: matriz diagonal com ganhos derivativos (amortecimento)
Kd = diag([30, 30, 30]);
% lambda: fator de amortecimento para cálculo da pseudo-inversa do Jacobiano (reduz efeitos numéricos ruins)
lambda = 0.1;

% Tempo de simulação. Define o passo de integração e o tempo total de simulação, criando um vetor de tempo discreto.
dt = 0.01;
T = 5;
tempo = 0:dt:T;

% --- Critérios de parada ---
posTol = 1e-3; % tolerância de posição (m)
velTol = 5e-3; % tolerância de velocidade cartesiana (m/s)
tMaxSegmento = 12; % tempo máximo por segmento (s) para evitar loop infinito

% Condições iniciais. Define as condições iniciais para as posições (ângulos) e velocidades das juntas.
theta = [0.5; -0.2; 0.2];
dtheta = [0.0; 0.0; 0.0];

% Posição desejada. Define o vetor de posição desejada do efetor final no espaço cartesiano.
Xd = [0.5; -0.2; 0.1];

% armazenar os thetas finais, para iniciar o próximo movimento dessa
% posição final

% --- Históricos globais (acumulam entre segmentos) ---
time_hist = []; % tempo global
t_total = 0; % relógio global

% Histórico para plot
% Inicializa matrizes para armazenar dados ao longo do tempo: posição do efetor final, ângulos das juntas e torque aplicado.
Xe_hist = zeros(3, length(tempo));
theta_hist = zeros(3, length(tempo));
torque_hist = zeros(3, length(tempo));

% Crie figura com dois subplots
figure('Color','w');
subplot(2,3,1); % Gráfico 3D do robô
ax1 = gca;
view(ax1, 0, 90); % Visão de cima
title('Visualização 3D (top view)');
axis equal;
hold on;

subplot(2,3,4); % Trajetória 2D
ax2 = gca;
title('Trajetória do efetuador');
xlabel('X (m)');
zlabel('Z (m)');
axis equal;
grid on;
hold on;

subplot(2,3,5); % Trajetória 2D
title('Trajetória do efetuador 3D');

% Subplot 3: Vista lateral (XZ)
subplot(2,3,2);
ax4 = gca;
title('Lateral (XZ)');
xlabel('X (m)');
zlabel('Z (m)');
axis equal;
grid on;
view(ax4, 0, 0); % Visão lateral pura
hold on;

axTorque = subplot(2,3,6);
hold(axTorque,'on');
h1 = plot(axTorque, tempo(1), torque_hist(1,1), 'r', 'LineWidth', 1.5);
h2 = plot(axTorque, tempo(1), torque_hist(2,1), 'g', 'LineWidth', 1.5);
h3 = plot(axTorque, tempo(1), torque_hist(3,1), 'b', 'LineWidth', 1.5);
xlabel(axTorque,'Tempo (s)');
ylabel(axTorque,'Torque (N·m)');
title(axTorque,'Torque em cada junta');
legend(axTorque, 'Junta 1','Junta 2','Junta 3');
grid(axTorque,'on');

% % Limites visuais iniciais
% xlim(ax1, [-1 1]); ylim(ax1, [-0.6 0.6]);
% xlim(ax3, [-1 1]); zlim(ax3, [0 1]);
% axis(ax2, [-1 1 -1 1 0 1]);

% --- Margens de segurança para checagem de espaço ---
eps_xy = 1e-3; % margem radial XY contra singularidade (m)
eps_z = 1e-3; % margem de Z contra topo/fundo (m)

% --- Limites derivados do robô ---
rmin = max(0, abs(l1 - l2) + eps_xy); % raio mínimo alcançável (XY)
rmax = (l1 + l2) - eps_xy; % raio máximo alcançável (XY)
zmin = (z0 - l3) + eps_z; % limite inferior de Z
zmax = (z0 - 0) - eps_z; % limite superior de Z (próximo ao topo)

fprintf('Espaço aproximado:\n');
fprintf(' XY: r ∈ [%.3f, %.3f] m\n', rmin, rmax);
fprintf(' Z : ∈ [%.3f, %.3f] m\n\n', zmin, zmax);

fprintf('Simulação interativa.\nEx.: [0.45 -0.15 0.12] | Digite q ou Enter para sair.\n\n');

%% --- Loop interativo ---
while true
    % ==== Entrada do usuário ====
    str = input('Digite o destino [X Y Z] em metros (ou q/Enter p/ sair): ','s');
    if isempty(str) || strcmpi(str,'q')
        fprintf('Encerrando.\n');
        break;
    end
    dest = str2num(str); %#ok<ST2NM>
    if numel(dest) ~= 3
        fprintf('Entrada inválida. Use três números: ex. 0.45 -0.15 0.12\n');
        continue;
    end
    Xd = dest(:);

    % ===== CHECAGEM XY =====
    r = hypot(Xd(1), Xd(2)); % distância radial ao eixo da base
    if r < rmin || r > rmax
        r_new = min(max(r, rmin), rmax); % "clip" radial
        if r > 0
            scale = r_new / r; % preserva direção (ângulo), ajusta módulo
            Xd(1:2) = Xd(1:2) * scale;
        else
            % caso especial: alvo exatamente na base (r = 0)
            Xd(1) = r_new; Xd(2) = 0;
        end
        warning('Destino XY fora do alcance. Ajustado para (%.3f, %.3f) (r=%.3f).', Xd(1), Xd(2), r_new);
    end

    % ===== CHECAGEM Z =====
    if Xd(3) < zmin || Xd(3) > zmax
        Xd(3) = min(max(Xd(3), zmin), zmax);
        warning('Destino Z fora do alcance. Ajustado para Z=%.3f (limites [%.3f, %.3f]).', Xd(3), zmin, zmax);
    end

    % --- Simulação deste segmento até atingir o alvo ou tempo máximo ---
    t_seg = 0;
    atingiu = false;

    k = 1; % inicializa o índice de contagem dentro do loop

    while t_seg < tMaxSegmento

%for k = 1:length(tempo)

    % Cinemática direta. Calcula a posição do efetor final no espaço cartesiano a partir dos ângulos das juntas 𝜃
    x = l1.*cos(theta(1)) + l2.*cos(theta(1)+theta(2));
    y = l1.*sin(theta(1)) + l2.*sin(theta(1)+theta(2));
    z = z0 - theta(3);
    X = [x; y; z];

    % Jacobiano. Calcula o Jacobiano geométrico do robô. O Jacobiano representa a derivada da posição do efetuador em relação às juntas.
    J11 = -l1.*sin(theta(1)) - l2.*sin(theta(1)+theta(2));
    J12 = -l2.*sin(theta(1)+theta(2));
    J13 = 0;
    J21 = l1.*cos(theta(1)) + l2.*cos(theta(1)+theta(2));
    J22 = l2.*cos(theta(1)+theta(2));
    J23 = 0;
    J31 = 0;
    J32 = 0;
    J33 = -1;
    J = [J11 J12 J13; J21 J22 J23; J31 J32 J33];

    % Erros
    % deltaX: erro de posição entre a posição desejada e a posição atual
    deltaX = Xd - X;
    % dX: velocidade atual do efetuador (obtida via Jacobiano)
    dX = J*dtheta;
    % delta_dX: erro de velocidade desejada, aqui é feito para frear o movimento atual (controle PD)
    delta_dX = -dX;

    % Controle PD no espaço operacional
    % Calcula a aceleração desejada no espaço operacional usando controle PD, que corrige erros de posição e velocidade.
    % Essa aceleração será convertida para velocidades articulares via Jacobiano.
     aX = Kp*deltaX + Kd*delta_dX;

    % Pseudo-inversa com amortecimento
    % Calcula a velocidade das juntas desejada para atingir aceleração aX, usando a pseudo-inversa do Jacobiano com amortecimento (damped least squares).
    % Amortecimento evita problemas de singularidade do Jacobiano e melhora estabilidade numérica.
    % Pseudo-inversa permite mapear velocidades do espaço operacional para velocidades articulares, mesmo para sistemas redundantes ou próximos de singularidades.
     J_damp = (J')*pinv(J*(J') + (lambda.^2).*eye(3));
     dtheta_d = J_damp*aX;

    % Dinâmica do robô RR
    % Termos de acoplamento (Coriolis)
    % Calcula a matriz de Coriolis e termos de acoplamento dinâmica.
    h = -m2.*l1.*l2.*sin(theta(2));
    C = [h.*dtheta(2), h.*sum(dtheta(1:2)), 0; -h.*dtheta(1), 0, 0; 0, 0, 0];

    % Matriz de inércia (3x3) que relaciona aceleração angular aos torques aplicadas, considerando inércia dos elos e efeito de configuração.
    M11 = I1 + I2 + m2.*l1.^2 + 2.*m2.*l1.*l2.*cos(theta(2));
    M12 = I2 + m2.*l1.*l2.*cos(theta(2));
    M21 = M12;
    M22 = I2;
    M = [M11 M12 0; M21 M22 0; 0 0 m3];
    % G é o vetor de forças gravitacionais atuando nas juntas (aqui só junta 3 tem peso vertical).
    G = [0; 0; m3*g];

    %Força virtual no espaço cartesiano (PD)
    %F = Kp*deltaX + Kd*delta_dX;

    %Torque via Jacobiano transposto
    %tau = J'*F;  

    Lambda = pinv(J*(M\J')); % robusto perto de singularidade
    aX_ref = Kp*deltaX + Kd*delta_dX;
    F = Lambda * aX_ref;
    tau = J'*F + C*dtheta + G; % compensa dinâmica em junta

    % Cálculo do torque e atualização das dinâmicas
    % Calcula o torque necessário para obter a aceleração desejada considerando dinâmica completa (inércia, Coriolis e gravidade).
     %tau = M*dtheta_d + C*dtheta + G;

    % Resolve a equação dinâmica para encontrar a aceleração angular real
    % resultante do torque aplicado. É a aceleração efetiva dada a torque e as forças dinâmicas.
    ddtheta = M\(tau - C*dtheta - G);

    % Integração
    % Atualiza as velocidades e posições das juntas pela integração numérica (método de Euler simples).
    dtheta = dtheta + ddtheta.*dt;
    theta = theta + dtheta.*dt;

    % Armazenar histórico
    Xe_hist(:, k) = X;
    theta_hist(:, k) = theta;
    torque_hist(:, k) = tau;

    robotConfig(1).JointPosition = theta(1);
    robotConfig(2).JointPosition = theta(2);
    robotConfig(3).JointPosition = theta(3);

    % Atualiza subplot 1: robô
    subplot(2,3,1);
    show(robot, robotConfig, 'PreservePlot', false, 'Frames','off','Parent',ax1);
    view(ax1, 0, 90); % Visão de cima
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('Vista superior (XY)');
    xlim(ax1, [-1 1]); ylim(ax1, [-0.6 0.6]);
    grid on;

    % Atualiza subplot 2: trajetória
    subplot(2,3,4);
    plot3(Xe_hist(1,:), Xe_hist(2,:), Xe_hist(3,:),'bo', 'MarkerSize',1, 'LineWidth', 1);
    plot3(Xd(1), Xd(2), Xd(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('Trajetória do efetuador');
    grid on; axis([-1 1 -1 1 0 1]); view(3);

    % Atualiza subplot 3: robô 3D
    subplot(2,3,5);
    show(robot, robotConfig, 'Visuals', 'on', 'Frames', 'off');
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('Visualização 3D');
    % Defina a visibilidade dos sistemas de referência como 'off'
    set(findall(gca, 'Tag', 'URDFLinkFrame'), 'Visible', 'off');
    axis([-1 1 -1 1 0 1]);view(3);

    % Subplot 4: lateral (XZ)
    subplot(2,3,2);
    show(robot, robotConfig, 'PreservePlot', false, 'Frames','off', 'Parent', ax4);
    view(ax4, 0, 0); % lateral
    title('Vista lateral (XZ)');
    xlim(ax4, [-1 1]); zlim(ax4, [0 1]);

    % Atualiza os dados do gráfico de torque
    if k > 1
        set(h1, 'XData', tempo(1:k), 'YData', torque_hist(1,1:k));
        set(h2, 'XData', tempo(1:k), 'YData', torque_hist(2,1:k));
        set(h3, 'XData', tempo(1:k), 'YData', torque_hist(3,1:k));
    end
    
    drawnow;
    %pause(0.0001);
%end

    % Atualiza tempo
    t_total = t_total + dt;
    t_seg = t_seg + dt;
    % Incrementa índice e tempo
    k = k + 1;

    % Critério de parada (alvo atingido)
        if norm(deltaX) < posTol && norm(dX) < velTol
            atingiu = true;
            break;
        end
    end

    if atingiu
        fprintf('Alvo (%.3f, %.3f, %.3f) atingido em %.2f s.\n Angulos (%.3f, %.3f, %.3f)', Xd(1),Xd(2),Xd(3), t_seg, theta(1)*180/pi,theta(2)*180/pi,theta(3)*180/pi);
    else
        fprintf('Não convergiu ao alvo em %.1f s (tente outro ponto ou ajuste ganhos).\n', tMaxSegmento);
    end
end

fprintf('Fim da simulação.\n');
% 
% % Gráfico de torque
% figure('Color','w');
% plot(tempo, torque_hist(1,:), 'r', 'LineWidth', 1.5); hold on;
% plot(tempo, torque_hist(2,:), 'g', 'LineWidth', 1.5);
% plot(tempo, torque_hist(3,:), 'b', 'LineWidth', 1.5);
% xlabel('Tempo (s)'); ylabel('Torque (N·m)');
% title('Torque em cada junta');
% legend('Junta 1','Junta 2','Junta 3');
% grid on;


% subplot(1,2,1);
% plot(Xe_hist(1,:), Xe_hist(2,:), 'b', 'LineWidth', 1.5); hold on;
% plot(Xd(1), Xd(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% xlabel('x [m]'); ylabel('y [m]');
% title('Trajetória no Espaço Cartesiano');
% legend('Trajetória do Efetuador', 'Objetivo');
% grid on; axis equal;
% subplot(1,2,2)
% plot(tempo,Xe_hist(1,:),tempo,Xe_hist(2,:));
% title('Comportamento Cartesiano no Tempo');
% legend('Posição: x(t)', 'Posição: y(t)');
% valor_min=min(min(Xe_hist(1,:)),min(Xe_hist(2,:)));
% valor_max=max(max(Xe_hist(1,:)),max(Xe_hist(2,:)));
% grid on; xlabel('Tempo [s]'); ylabel('Amplitude [m]');
% axis([0 tempo(1,end) valor_min valor_max]);
% set(gcf,'Color',[1 1 1]);

% % Plot final
% figure;
% plot3(Xe_hist(1,:), Xe_hist(2,:), Xe_hist(3,:), 'b', 'LineWidth', 2); hold on;
% plot3(Xd(1), Xd(2), Xd(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
% title('Trajetória do efetuador com dinâmica completa nos 3 elos');
% grid on; axis equal; view(3);
