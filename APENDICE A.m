robot = importrobot('omronEcobra600.urdf');
robotConfig = homeConfiguration(robot);

l1 = 0.325; 
l2 = 0.275; 
l3 = 0.210; 

z0 = 0.387;

m1 = 2.0; 
m2 = 2.0; 
m3 = 0.4;

I1 = 0.01; 
I2 = 0.01; 
I3 = 0.001;

g = 9.81;

Kp = diag([120, 120, 120]);
Kd = diag([30, 30, 30]);
lambda = 0.1;

dt = 0.01;
T = 5;
tempo = 0:dt:T;

posTol = 1e-3;
velTol = 5e-3;
tMaxSegmento = 12;

theta = [0.5; -0.2; 0.2];
dtheta = [0.0; 0.0; 0.0];

t_total = 0;     

Xe_hist = zeros(3, length(tempo));
theta_hist = zeros(3, length(tempo));
torque_hist = zeros(3, length(tempo));

figure('Color','w');
subplot(2,3,1);
ax1 = gca;
view(ax1, 0, 90);
title('Visualização 3D (top view)');
axis equal;
hold on;

subplot(2,3,4);
ax2 = gca;
title('Trajetória do efetuador');
xlabel('X (m)');
zlabel('Z (m)');
axis equal;
grid on;
hold on;

subplot(2,3,5);
title('Trajetória do efetuador 3D');

subplot(2,3,2);
ax4 = gca;
title('Lateral (XZ)');
xlabel('X (m)');
zlabel('Z (m)');
axis equal;
grid on;
view(ax4, 0, 0);
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

eps_xy = 1e-3;
eps_z  = 1e-3;

rmin = max(0, abs(l1 - l2) + eps_xy);
rmax = (l1 + l2) - eps_xy;           
zmin = (z0 - l3) + eps_z;              
zmax = (z0 - 0)  - eps_z;              

fprintf('Espaço aproximado:\n');
fprintf('  XY: r ∈ [%.3f, %.3f] m\n', rmin, rmax);
fprintf('  Z :   ∈ [%.3f, %.3f] m\n\n', zmin, zmax);

fprintf(['Simulação interativa.\nEx.: [0.45 -0.15 0.12]' ...
    '  |  Digite q ou Enter para sair.\n\n']);

%% --- Loop interativo ---
while true
    str = input(['Digite o destino [X Y Z] ' ...
        'em metros (ou q/Enter p/ sair): '],'s');
    if isempty(str) || strcmpi(str,'q')
        fprintf('Encerrando.\n');
        break;
    end
    dest = str2num(str); %#ok<ST2NM>
    if numel(dest) ~= 3
        fprintf(['Entrada inválida. ' ...
            'Use três números: ex. 0.45 -0.15 0.12\n']);
        continue;
    end

    Xd = dest(:);

    r = hypot(Xd(1), Xd(2));
    if r < rmin || r > rmax
        r_new = min(max(r, rmin), rmax);
        if r > 0
            scale = r_new / r;
            Xd(1:2) = Xd(1:2) * scale;
        else
            Xd(1) = r_new; Xd(2) = 0;
        end
        warning(['Destino XY fora do alcance. ' ...
            'Ajustado para (%.3f, %.3f) (r=%.3f).' ...
            ], Xd(1), Xd(2), r_new);
    end

    if Xd(3) < zmin || Xd(3) > zmax
        Xd(3) = min(max(Xd(3), zmin), zmax);
        warning(['Destino Z fora do alcance. ' ...
            'Ajustado para Z=%.3f (limites ' ...
            '[%.3f, %.3f]).'], Xd(3), zmin, zmax);
    end

    t_seg = 0;
    atingiu = false;

    k = 1;

    while t_seg < tMaxSegmento
        
        x = l1.*cos(theta(1)) + l2.*cos(theta(1)+theta(2));
        y = l1.*sin(theta(1)) + l2.*sin(theta(1)+theta(2));
        z = z0 - theta(3);
        X = [x; y; z];
    
        J11 = -l1.*sin(theta(1)) - l2.*sin(theta(1)+theta(2));
        J12 = -l2.*sin(theta(1)+theta(2));
        J13 =  0;
        J21 =  l1.*cos(theta(1)) + l2.*cos(theta(1)+theta(2));
        J22 =  l2.*cos(theta(1)+theta(2));
        J23 =  0;
        J31 =  0;
        J32 =  0;
        J33 = -1;
        J = [J11 J12 J13; J21 J22 J23; J31 J32 J33];
    
        deltaX = Xd - X;
        dX = J*dtheta;
        delta_dX = -dX;

        aX = Kp*deltaX + Kd*delta_dX;
    
        J_damp = (J')*pinv(J*(J') + (lambda.^2).*eye(3));
        dtheta_d = J_damp*aX;
    
        h = -m2.*l1.*l2.*sin(theta(2));
        C = [h.*dtheta(2), h.*sum(dtheta(1:2)), 0; ...
            -h.*dtheta(1), 0, 0; 0, 0, 0];
    
        M11 = I1 + I2 + m2.*l1.^2 + 2.*m2.*l1.*l2.*cos(theta(2));
        M12 = I2 + m2.*l1.*l2.*cos(theta(2));
        M21 = M12;
        M22 = I2;
        M = [M11 M12 0; M21 M22 0; 0 0 m3];
        G = [0; 0; m3*g];           
    
        Lambda  = pinv(J*(M\J'));
        aX_ref  = Kp*deltaX + Kd*delta_dX;
        F       = Lambda * aX_ref;
        tau     = J'*F + C*dtheta + G;

        ddtheta = M\(tau - C*dtheta - G);
    
        dtheta = dtheta + ddtheta.*dt;
        theta = theta + dtheta.*dt;

        Xe_hist(:, k) = X;
        theta_hist(:, k) = theta;
        torque_hist(:, k) = tau;
    
        robotConfig(1).JointPosition = theta(1);
        robotConfig(2).JointPosition = theta(2);
        robotConfig(3).JointPosition = theta(3);
    
        subplot(2,3,1);
        show(robot, robotConfig, 'PreservePlot', ...
            false, 'Frames','off','Parent',ax1);
        view(ax1, 0, 90); 
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title('Vista superior (XY)');
        xlim(ax1, [-1 1]); ylim(ax1, [-0.6 0.6]);
        grid on;
    
        subplot(2,3,4);
        plot3(Xe_hist(1,:), Xe_hist(2,:), Xe_hist(3,:), ...
            'bo', 'MarkerSize',1, 'LineWidth', 1);
        plot3(Xd(1), Xd(2), Xd(3), 'ro', 'MarkerSize' ...
            , 10, 'LineWidth', 2);
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title('Trajetória do efetuador');
        grid on; axis([-1 1 -1 1 0 1]); view(3);
    
        subplot(2,3,5);
        show(robot, robotConfig, 'Visuals', 'on', ...
            'Frames', 'off');
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title('Visualização 3D');
        set(findall(gca, 'Tag', 'URDFLinkFrame'), ...
            'Visible', 'off');
        axis([-1 1 -1 1 0 1]);view(3);
    
        subplot(2,3,2);
        show(robot, robotConfig, 'PreservePlot', false, ...
            'Frames','off', 'Parent', ax4);
        view(ax4, 0, 0);
        title('Vista lateral (XZ)');
        xlim(ax4, [-1 1]); zlim(ax4, [0 1]);
    
        if k > 1
            set(h1, 'XData', tempo(1:k), 'YData', ...
                torque_hist(1,1:k));
            set(h2, 'XData', tempo(1:k), 'YData', ...
                torque_hist(2,1:k));
            set(h3, 'XData', tempo(1:k), 'YData', ...
                torque_hist(3,1:k));
        end
        
        drawnow;
    
        t_total = t_total + dt;
        t_seg   = t_seg   + dt;

        k = k + 1;
    
        if norm(deltaX) < posTol && norm(dX) < velTol
            atingiu = true;
            break;
        end
    end

    if atingiu
        fprintf(['Alvo (%.3f, %.3f, %.3f) atingido ' ...
            'em %.2f s.\n Angulos (%.3f, %.3f, %.3f)' ...
            ], Xd(1),Xd(2),Xd(3), t_seg, ...
            theta(1)*180/pi,theta(2)*180/pi, ...
            theta(3)*180/pi);
    else
        fprintf(['Não convergiu ao alvo em %.1f s ' ...
            '(tente outro ponto ou ajuste ganhos).' ...
            '\n'], tMaxSegmento);
    end
end

fprintf('Fim da simulação.\n');