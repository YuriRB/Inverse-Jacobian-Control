%% Trajetória 
%% Alunos: Reginaldo Barbosa/Yuri Ruzzene
      

%%Matrizes de Posição e Rotação do Robô Scara Epson T3 401SS

clear,close all, clc

%%syms cria icógnitas para gerar as funções 
%%D3z = distancia de p3 em relação a p2 no eixo z
syms D3z teta1 teta2 teta4 real;

%% Matrizes de posição dos referenciais
P01 = [ 0; 0; 0.2098];
P12 = [ 0; 0.225; 0.035];
P23 = [ 0; 0.175; D3z];
P34 = [0; 0 ; -0.093];

%% Matrizes de rotação
R01 = [cos(teta1) -sin(teta1) 0; 
       sin(teta1)  cos(teta1) 0; 
       0               0      1];
   

R12 = [cos(teta2) -sin(teta2) 0; 
        sin(teta2) cos(teta2) 0; 
        0             0       1];


R23 = [1 0 0; 
       0 1 0; 
       0 0 1];


R34 = [cos(teta4) -sin(teta4) 0; 
      -sin(teta4) -cos(teta4) 0; 
            0           0     -1];



%% Matrizes T
T01 = [R01 P01; 0 0 0 1];
T12 = [R12 P12; 0 0 0 1];
T23 = [R23 P23; 0 0 0 1];
T34 = [R34 P34; 0 0 0 1];

%% Matriz T do ponto final do robô em relação ao seu referencial fixo
T04 = (T01 * (T12 * (T23 * T34))) ;

nT04 = subs(T04, [teta1, teta2, teta4, D3z], [ deg2rad(180) deg2rad(90) 0 0]);


%%  Lógica inversa Método 1
 syms px py pz real;
% 
% eq_x = T04(1,4) == px;
% eq_y = T04(2,4) == py;
% eq_z = T04(3,4) == pz;
% 
% t1 = isolate(eq_x,teta1)
% t2 = isolate(eq_y,teta2)
% p3 = isolate(eq_z,D3z)

%% Método NR 
F = real(T04(1:3,4) - [px py pz]');
dF = jacobian(F,[teta1 teta2 D3z]);

teta = [0 0 0]';%P=[0 400 151.8]' P2=[175 -225 151.8]

t1 = teta(1,1); t2=teta(2,1); D3 = teta(3,1);

p = [175 225 100]';%ponto onde tudo é 0

x = p(1,1); y=p(2,1);z=p(3,1);

Fp=[1 1 1];
tolerancia=1e-9;
iteracao=0;

while tolerancia<=max(abs(Fp));
 
 Fp=[-175*cos(t1)*sin(t2)-sin(t1)*(175*cos(t2)+225)-x;
     cos(t1)*(175*cos(t2)+225)-175*sin(t1)*sin(t2)-y;
     D3+151.8-z];

dFp=[175*sin(t1)*sin(t2)-cos(t1)*(175*cos(t2)+225),175*sin(t1)*sin(t2)-175*cos(t1)*cos(t2),0;
-175*cos(t1)*sin(t2)-sin(t1)*(175*cos(t2)+225),-175*cos(t1)*sin(t2)-175*cos(t2)*sin(t1),0;
0,0,1];

teta=teta-pinv(dFp)*Fp;
t1=teta(1,1);
t2=teta(2,1);
D3=teta(3,1);
iteracao=iteracao+1;

if iteracao > 100
    break;
end

end


Fpo=[-175.*cos(t1).*sin(t2)-sin(t1).*(175.*cos(t2)+225);
     cos(t1).*(175.*cos(t2)+225)-175.*sin(t1).*sin(t2);
     D3+(759./5)]

format short
Resultado = [rad2deg(t1)-(floor(rad2deg(t1)/360)*360) rad2deg(t2)-(floor(rad2deg(t2)/360)*360) D3]'

%% Trajetoria
close all
a0 = [];
a1 = [];
a2 = [];
a3 = [];

theta0 = [];
thetaF = [];

tf = 20;
T = 0.1;

%Junta 1
theta0(1) = [0];
thetaF(1) = [pi];

a0(1) = theta0(1);
a1(1) = 0;
a2(1) = 3*(thetaF(1) - theta0(1))/(tf^2);
a3(1) = -2*(thetaF(1) - theta0(1))/(tf^3);

%Junta 2

theta0(2) = [0];
thetaF(2) = [-pi];

a0(2) = theta0(2);
a1(2) = 0;
a2(2) = 3*(thetaF(2) - theta0(2))/(tf^2);
a3(2) = -2*(thetaF(2) - theta0(2))/(tf^3);

%junta 3

theta0(3) = [0];
thetaF(3) = [0.15];

a0(3) = theta0(3);
a1(3) = 0;
a2(3) = 3*(thetaF(3) - theta0(3))/(tf^2);
a3(3) = -2*(thetaF(3) - theta0(3))/(tf^3);

for i = 1: tf/T + 1
    t = (i-1)*T;
    t_p(:,i)=t;
    
    %Junta 1
    d = a0(1) + a1(1).*t + a2(1).*t^2 + a3(1).*t^3;
    dd= a1(1) + 2.*a2(1).*t + 3.*a3(1).*t^2;
    ddd= 2.*a2(1) + 6.*a3(1).*t;
    
    d_p1(:,i)=d;
    v_p1(:,i)=dd;
    a_p1(:,i)=ddd;
    
    %Junta 2
    d = a0(2) + a1(2).*t + a2(2).*t^2 + a3(2).*t^3;
    dd= a1(2) + 2.*a2(2).*t + 3.*a3(2).*t^2;
    ddd= 2.*a2(2) + 6.*a3(2).*t;
    
    d_p2(:,i)=d;
    v_p2(:,i)=dd;
    a_p2(:,i)=ddd;
    
    %Junta 3
    d = a0(3) + a1(3).*t + a2(3).*t^2 + a3(3).*t^3;
    dd= a1(3) + 2.*a2(3).*t + 3.*a3(3).*t^2;
    ddd= 2.*a2(3) + 6.*a3(3).*t;
    
    v_p3(:,i)=dd;
    a_p3(:,i)=ddd;    d_p3(:,i)=d;

    
    
    %Cinemática direta
    
    nT04 = subs(T04, [teta1, teta2, D3z, teta4], [ d_p1(:,i)  d_p2(:,i) d_p3(:,i) 0]);
    
    plot3(nT04(1,4),nT04(2,4),nT04(3,4),'o'),axis([-420 410 -420 420 -400 400]), xlabel('x'), ylabel('y'), zlabel('z');
    hold on;
    
end

% 
% figure
% subplot(1,3,1), plot(d_p1);
% subplot(1,3,2), plot(v_p1);
% subplot(1,3,3), plot(a_p1);
% 
% figure
% subplot(1,3,1), plot(d_p2);
% subplot(1,3,2), plot(v_p2);
% subplot(1,3,3), plot(a_p2);
% 
% figure
% subplot(1,3,1), plot(d_p3);
% subplot(1,3,2), plot(v_p3);
% subplot(1,3,3), plot(a_p3);
% 
% 
figure
subplot(1,3,1)
plot(t_p,d_p1(1,:),'LineWidth',1.2);
grid on;
xlabel("Tempo [s]")
ylabel("Posição [rad]")
subplot(1,3,2)
plot(t_p,v_p1,'LineWidth',1.2);
grid on;
title("Elo 1")
xlabel("Tempo [s]")
ylabel("Velocidade [rad/s]")
subplot(1,3,3)
plot(t_p,a_p1,'LineWidth',1.2);
grid on;
xlabel("Tempo [s]")
ylabel("Aceleração [rad/s^2]")


figure
subplot(1,3,1)
plot(t_p,d_p2,'LineWidth',1.2);
grid on;
xlabel("Tempo [s]")
ylabel("Posição [rad]")
subplot(1,3,2)
plot(t_p,v_p2,'LineWidth',1.2);
grid on;
title("Elo 2")
xlabel("Tempo [s]")
ylabel("Velocidade [rad/s]")
subplot(1,3,3)
plot(t_p,a_p2,'LineWidth',1.2);
grid on;
xlabel("Tempo [s]")
ylabel("Aceleração [rad/s^2]")

figure
subplot(1,3,1)
plot(t_p,d_p3,'LineWidth',1.2);
grid on;
xlabel("Tempo [s]")
ylabel("Posição [rad]")
subplot(1,3,2)
plot(t_p,v_p3,'LineWidth',1.2);
grid on;
title("Elo 3")
xlabel("Tempo [s]")
ylabel("Velocidade [rad/s]")
subplot(1,3,3)
plot(t_p,a_p3,'LineWidth',1.2);
grid on;
xlabel("Tempo [s]")
ylabel("Aceleração [rad/s^2]")

%% Dinâmica
% 
% %peso do robô total 16,5kg
% syms g real
% P1_c1 = [0 0.1125 0]';
% R10 = R01';
% dv00 = [ 0 0 -g]';
% w00 = [0 0 0]';
% dw00 = [0 0 0]';
% Z=[0 0 1]';
% m1 = 1.2;
% 
% 
% % W(i+1)(i+1) = R(i+1)i . Wii + theta_ponto_i+1 . z(i+1)(i+1);
% % dW(i+1)(i+1) = R(i+1)i . Wii + R(i+1)i . Wii x theta_ponto_i+1 . z(i+1)(i+1) + theta_2pontos(i+1) . z(i+1)(i+1); 
% % dv(i+1)(i+1) = R(i+1)i . (dWii x Pi(i+1) + wii x (wii x pi(i+1) + dVii))
% % dVc(i+1))(i+1) = dW(i+1)(i+1) x Pc(i+1)(i+1) + w(i+1)(i+1) x (w(i+1)(i+1)xPc(i+1)(i+1) + dv(i+1)(i+1)
% % F(i+1)(i+1) = m(i+1) . dVc(i+1)(i+1)
% % N(i+1)(i+1) = Ic(i+1)(i+1). dW(i+1)(i+1) + w(i+1)(i+1) xIc(i+1)(i+1).w(i+1)(i+1)
% 
% % x= k
% % y= l
% % z= h 
% 
% 
% %i = 0
% %x = 0.090 y = 0.225 z = 0.035
% Ic11 = [ (0.035^2 + 0.225^2)*m1/12             0                             0;
%                       0             (0.035^2 + 0.090^2)*m1/12                0;
%                       0                        0               (0.090^2 + 0.225^2)*m1/12 ];
% 
% syms dhteta1 ddhteta1 real 
% 
% w11 = R10 * w00 + dhteta1 * Z;
% dw11 = R10 * dw00 + cross(R10 * w00, dhteta1 * Z) + ddhteta1 * Z; 
% dv11 = R10 * (cross(dw00,P01) + cross(w00,cross(w00, P01)) + dv00);
% dvc11 = cross(dw11,P1_c1) + cross(w11,cross(w11,P1_c1)) + dv11;
% F11 = m1 * dvc11;
% N11 = Ic11 * dw11 + cross(w11,Ic11 * w11);
% 
% 
% 
% 
% % x= k
% % y= l
% % z= h 
% 
% %x = 0.129 y = 0.175 z = 0.1765
% %i = 1
% P2_c2 = [ 0 0.0875 0]';
% m2 = 4;
% syms dhteta2 ddhteta2 real 
% Ic22 = [ (0.175^2 + 0.1765^2)*m2/12             0                             0;
%                       0             (0.1295^2 + 0.1765^2)*m2/12                0;
%                       0                        0               (0.175^2 + 0.1295^2)*m2/12 ];
% w22 = R10 * w11 + dhteta2 * Z;
% dw22 = R10 * dw11 + cross(R10 * w11, dhteta2 * Z) + ddhteta2 * Z; 
% dv22 = R10 * (cross(dw11,P12) + cross(w11,cross(w11, P12)) + dv11);
% dvc22 = cross(dw22,P2_c2) + cross(w22,cross(w22,P2_c2)) + dv22;
% F22 = m2 * dvc22;
% N22 = Ic22 * dw22 + cross(w22,Ic22 * w22);
% 
% % x= k
% % y= l
% % z= h 
% 
% %x = 0.00857 y = 0.00857 z = 0.3407
% %i = 2
% P3_c3 = [ 0 0 0.17035]';
% m3 = 0.4;
% syms dhteta3 ddhteta3 real
% 
% Ic33 = [ (0.00857^2 + 0.3407^2)*m3/12             0                         0;
%                       0                  (0.00857^2 + 0.3407^2)*m3/12       0;
%                       0                             0                    (0.00857^2 + 0.00857^2)*m3/12   ];
% 
% w33 = R10 * w22 + dhteta3 * Z;
% dw33 = R10 * dw22 + cross(R10 * w22, dhteta3 * Z) + ddhteta3 * Z; 
% dv33 = R10 * (cross(dw22,P23) + cross(w22,cross(w22, P23)) + dv22);
% dvc33 = cross(dw22,P3_c3) + cross(w33,cross(w33,P3_c3)) + dv33;
% F33 = m3 * dvc33;
% N33 = Ic33 * dw33 + cross(w33,Ic33 * w33);
% 
% %RETORNO
% % i = 3
% 
% f44 = [0 0 0]';
% n44 = [0 0 0]';
% tau4 = 0;
% 
% f33 = R34'*f44+F33;
% n33 = N33+R34'*n44+cross(P3_c3,F33)+cross(P34,R34'*f44);
% tau3 = f33'*Z;
% 
% % i = 2
% f22 = R23'*f33+F22;
% n22 = N22+R23'*n33+cross(P2_c2,F22)+cross(P23,R23'*f33);
% tau2 = n22'*Z;
% 
% % i = 1
% f11 = R12'*f22+F11;
% n11 = N11+R12'*n22+cross(P1_c1,F11)+cross(P23,R23'*f22);
% tau1 = n11'*Z;
% 
% %Taxa de amostragem/tempo de simulação
% Tx=1e-4; seg=10;
% t1=0; t2=0; t3=0; dt1=0; dt2=0; dt3=0;
% 
% for i=1:seg/Tx+1
%        t=Tx.*(i-1);
%        t_p(:,i)=t;
%        
%        ta1= heaviside(t)-heaviside(t-0.1);
%        ta2= (heaviside(t)-heaviside(t-0.1));
%        ta3= (heaviside(t)-heaviside(t-0.1));
%        
%        Tau=[ta1; ta2; ta3]; Tau_p(:,i)=Tau;
%       
%        T=[t1; t2; t3];
%        T_p(:,i)=T;
%        
%        dT=[dt1; dt2; dt3];
%        dT_p(:,i)=dT;
%        
%        
%        M=[ 0.1514+0.2362*cos(t1)-0.0315*(sin(t1))^2+0.0245*cos(t1)+0.0315*(cos(t1))^2, 0.1252+0.0245*cos(t1), 1.4689e-5;
%            0.0639+0.0788*cos(t1)-0.0158*sin(t1)*sin(t1)+0.0123*cos(t1)+0.0158*cos(t1)^2,0.0639+0.0123*cos(t1),1.4689e-5; 
%            0,0,0];
%    
%        
%        V = [0.2362*sin(t1)*dt1^2+(7*sin(t1)*((dt1+1.4689*10^(-5)*dt2)*((7*dt1)/40+(7*dt2)/40)+(9*dt1^2*cos(t1))/40))/50+(63*dt1^2*cos(t1)*sin(t1))/2000;
%             0.0788*dt1^2*sin(t1)+0.0123*sin(t1)*dt1^2+0.0245*sin(t1)*dt1*dt2+0.0123*sin(t1)*dt2^2+0.0158*sin(t1)*dt1^2*cos(t1)+0.0158*cos(t1)*sin(t1)*dt1^2;
%             0] ;
%                                                                                   
%                                                                                   
%        G=[0;0;-173/25];
%        
%        F=2*[dt1; 
%                dt2;
%                dt3]; %Vetor de amortecimento
%        
%        ddT=pinv(M)*(Tau-V-G-F);
%        
%        dT=dT+Tx.*ddT;
%        
%        T=T+Tx.*dT;
%        t1=T(1,1);
%        t2=T(2,1);
%        t3=T(3,1);
%        dt1=dT(1,1);
%        dt2=dT(2,1);
%        dt3=dT(3,1);
%     
% end
% 
% 
% subplot(1,3,1), plot(t_p(1,:),T_p(1,:)), grid on
% xlabel('Tempo'), ylabel('\theta_1')
% subplot(1,3,2), plot(t_p(1,:),T_p(2,:)), grid on
% xlabel('Tempo'), ylabel('\theta_2')
% subplot(1,3,3), plot(t_p(1,:),T_p(3,:)), grid on
% xlabel('Tempo'), ylabel('\theta_3')
% set(gcf,'Color',[1 1 1])
% close all;
% %% CONTROLE NÃO LINEAR
% 
% %Aula 11 - Sistema com amortecedor não linear - Controle Não Linear - Ex.4
% 
% kp=1;
% kv=2;
% T=1e-3;
% t_max=10;
% 
% THETA1=zeros(3,1);
% THETA2=zeros(3,1);
% THETA3=zeros(3,1);
% 
% THETA1d=[1, 0, 0]'; %[Posição desejada, Velocidade desejada, Aceleração desejada
% THETA2d=[1, 0, 0]';
% THETA3d=[1, 0, 0]';
% 
% ThetaFinal = zeros(3,1);
% dThetaFinal = zeros(3,1);
% ddThetaFinal = zeros(3,1);
% 
% for i=1:t_max./T+1;
%     t=(i-1).*T;
%     t_p(1,i)=t;
%     
%     ep1 = THETA1d(1)-THETA1(1);
%     ev1 = THETA1d(2)-THETA1(2);
%     fl1 = THETA1d(3)+kv*ev1+kp*ep1;
%     
%     ep2 = THETA2d(1)-THETA2(1);
%     ev2 = THETA2d(2)-THETA2(2);
%     fl2 = THETA2d(3)+kv*ev2+kp*ep2;
%     
%     ep3 = THETA3d(1)-THETA3(1);
%     ev3 = THETA3d(2)-THETA3(2);
%     fl3 = THETA3d(3)+kv*ev3+kp*ep3;
%     
%     fl = [fl1;fl2;fl3];
%     
%     f = M*fl+V+G;
%     
%     %Atualizar M, V, e G
%     
%     M = [ 0.1514+0.2362*cos(THETA1(1))-0.0315*(sin(THETA1(1)))^2+0.0245*cos(THETA1(1))+0.0315*(cos(THETA1(1)))^2, 0.1252+0.0245*cos(THETA1(1)), 1.4689e-5;
%            0.0639+0.0788*cos(THETA1(1))-0.0158*sin(THETA1(1))*sin(THETA1(1))+0.0123*cos(THETA1(1))+0.0158*cos(THETA1(1))^2,0.0639+0.0123*cos(THETA1(1)),1.4689e-5; 
%            0,0,0];
%        
%     V = [0.2362*sin(THETA1(1))*THETA1(2)^2+(7*sin(THETA1(1))*((THETA1(2)+1.4689*10^(-5)*THETA2(2))*((7*THETA1(2))/40+(7*THETA2(2))/40)+(9*THETA1(2)^2*cos(THETA1(1)))/40))/50+(63*THETA1(2)^2*cos(THETA1(1))*sin(THETA1(1)))/2000;
%             0.0788*THETA1(2)^2*sin(THETA1(1))+0.0123*sin(THETA1(1))*THETA1(2)^2+0.0245*sin(THETA1(1))*THETA1(2)*THETA2(2)+0.0123*sin(THETA1(1))*THETA2(2)^2+0.0158*sin(THETA1(1))*THETA1(2)^2*cos(THETA1(1))+0.0158*cos(THETA1(1))*sin(THETA1(1))*THETA1(2)^2;
%             0] ;   
%        
%     G=[0;0;-173/25];
%     
%     ddThetaFinal = pinv(M)*(-V-G+f);
%     dThetaFinal = dThetaFinal+ddThetaFinal*T;
%     ThetaFinal = ThetaFinal+dThetaFinal*T;
%     
%     THETA1 = [ThetaFinal(1,1); dThetaFinal(1,1); ddThetaFinal(1,1)];
%     THETA2 = [ThetaFinal(2,1); dThetaFinal(2,1); ddThetaFinal(2,1)];
%     THETA3 = [ThetaFinal(3,1); dThetaFinal(3,1); ddThetaFinal(3,1)];
%     
%     THETA_p1(:,i) = ThetaFinal(1,1);
%     THETA_p2(:,i) = ThetaFinal(2,1);
%     THETA_p3(:,i) = ThetaFinal(3,1);
%     
%     dTHETA_p1(:,i) = dThetaFinal(1,1);
%     dTHETA_p2(:,i) = dThetaFinal(2,1);
%     dTHETA_p3(:,i) = dThetaFinal(3,1);
%     
%     ddTHETA_p1(:,i) = ddThetaFinal(1,1);
%     ddTHETA_p2(:,i) = ddThetaFinal(2,1);
%     ddTHETA_p3(:,i) = ddThetaFinal(3,1);
%     
%     
%     
% %     X_p(:,i)=X(1:2,:);
% %     x3_p(1,i)=x3;
% %     X
% %     x1=X_p(1,i); % posição
% %     x2=X_p(2,i);
% %     Xd_p(:,i)=Xd;
% %     x1d=Xd_p(1,i);
% %     x2d=Xd_p(2,i);
% %     x3d=Xd_p(3,i);
% %                                                                                               
% %     ep=x1d-x1;
% %     ev=x2d-x2;
% %     
% %     alpha=M(1,:);
% %     %beta=bc.*sign(x2)+q.*x1;
% %     beta=V(1,:)+G(1,:);
% %     
% %     fl=x3d+kp.*ep+kv.*ev;
% %     f=alpha.*fl+beta;
% %     %dX=[x2;-G/m-V/m]+[0;1/m].*f;
% %     %dX=[x2;-(q./m).*x1-(bc./m).*sign(x2)]+[0;1./m].*f;
% %     dX=pinv(M(:,1)).*(Tau-V(1,1)-G(1,1)-f);
% %     x3=dX(2,1);
% %     X=X(:,1)+dX(1:2,1).*T;
% end
% % 
% % subplot(1,3,1), hold on, plot(t_p(1,1:10001),X_p(1,:)), grid on
% % xlabel('Tempo, s')
% % ylabel('Posição')
% % subplot(1,3,2), hold on, plot(t_p(1,1:10001),X_p(2,:)), grid on
% % title('Controle Não Linear')
% % xlabel('Tempo, s')
% % ylabel('Velocidade')
% % subplot(1,3,3), hold on, plot(t_p(1,1:10001),x3_p(1,:)), grid on
% % xlabel('Tempo, s')
% % ylabel('Aceleração')
% % set(gcf,'Color',[1 1 1])
% % 
% % subplot(1,3,1), hold on, plot(t_p(1,1:10001),X_p(1,:)), grid on
% % xlabel('Tempo, s')
% % ylabel('Posição')
% % subplot(1,3,2), hold on, plot(t_p(1,1:10001),X_p(2,:)), grid on
% % title('Controle Não Linear')
% % xlabel('Tempo, s')
% % ylabel('Velocidade')
% % subplot(1,3,3), hold on, plot(t_p(1,1:10001),x3_p(1,:)), grid on
% % xlabel('Tempo, s')
% % ylabel('Aceleração')
% % set(gcf,'Color',[1 1 1])
% figure(1)
% subplot(1,3,1), hold on, plot(t_p(1,1:10001),THETA_p1(1,:)), grid on
% title('Posição Elo 1')
% xlabel('Tempo [s]')
% ylabel('Posição [rad]')
% subplot(1,3,2), hold on, plot(t_p(1,1:10001),THETA_p2(1,:)), grid on
% title('Posição Elo 2')
% xlabel('Tempo [s]')
% ylabel('Posição [rad]')
% subplot(1,3,3), hold on, plot(t_p(1,1:10001),THETA_p3(1,:)), grid on
% title('Posição Elo 3')
% xlabel('Tempo [s]')
% ylabel('Posição [rad]')
% 
% figure(2)
% subplot(1,3,1), hold on, plot(t_p(1,1:10001),dTHETA_p1(1,:)), grid on
% title('Velocidade Elo 1')
% xlabel('Tempo [s]')
% ylabel('Velocidade [rad/s]')
% subplot(1,3,2), hold on, plot(t_p(1,1:10001),dTHETA_p2(1,:)), grid on
% title('Velocidade Elo 2')
% xlabel('Tempo [s]')
% ylabel('Velocidade [rad/s]')
% subplot(1,3,3), hold on, plot(t_p(1,1:10001),dTHETA_p3(1,:)), grid on
% title('Velocidade Elo 3')
% xlabel('Tempo [s]')
% ylabel('Velocidade [rad/s]')
% figure(3)
% subplot(1,3,1), hold on, plot(t_p(1,1:10001),ddTHETA_p1(1,:)), grid on
% title('Aceleração Elo 1')
% xlabel('Tempo [s]')
% ylabel('Aceleração [rad/s^2]')
% subplot(1,3,2), hold on, plot(t_p(1,1:10001),ddTHETA_p2(1,:)), grid on
% title('Aceleração Elo 2')
% xlabel('Tempo [s]')
% ylabel('Aceleração [rad/s^2]')
% subplot(1,3,3), hold on, plot(t_p(1,1:10001),ddTHETA_p3(1,:)), grid on
% title('Aceleração Elo 3')
% xlabel('Tempo [s]')
% ylabel('Aceleração [rad/s^2]')

