clear all;
clc;

T_s = 0.01;
t = 0:T_s:25;
x = zeros(2,length(t));
x_dot = zeros(2,length(t));
y = zeros(1,length(t));
y(1) = x(1,1);

omega_n = 10;
zeta = 0.17;
K = 2;
s = tf('s');
G = K*omega_n^2/(s^2 + 2*zeta*omega_n*s + omega_n^2);   % Stable
% G = 2/(s^2 - 0.25*s + 12);                              % Unstable

omega = linspace(0,13.42,length(t));
u = sin(omega.*t);
[n,d] = tfdata(c2d(G,T_s),'v');
[A,B,C,D] = tf2ss(n,d);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultLineLineWidth', 2);

x_limit = 24.51;
u_offset = 15;
for i = 2:length(t)
    x(:,i) = A*x(:,i-1) + B*u(i-1);
    y(:,i) = C*x(:,i) + D*u(i);
    if (mod(i,10) == 0)
        figure(1);
        u_plot = u(i) + u_offset;
        scatter(0,u_plot,254,'sq','markerFaceColor','k','MarkerEdgeColor','k');
        hold on;
        plot([0 4],[u_plot, u_plot],'--k');
        [x_s,y_s] = spring(0,u_plot,0,y(i),10,154,0.01);
        plot(x_s,y_s,'LineWidth',2);
        plot(-t(i:-1:1),u(1:i)+u_offset,'k');
        plot(-t(i:-1:1),y(1:i),'b');
        scatter(0,y(i),254,'o','markerFaceColor','b','MarkerEdgeColor','b');
        plot([0 4],[y(i), y(i)],'--k');
        plot([2.5 2.5],[u_plot, y(i)],'k');
        text(-21,-11,strcat(['$\omega = $',num2str(omega(i)),' rad/s']),'fontsize',20);
        hold off;
        grid on;
        axis([-25 5.4 -21 21]);
        title('Input in black, output in blue');
    end
end
[up,lo] = envelope(y);
[up_u,lo_u] = envelope(u);
%%
figure(4);
subplot(2,2,[1 3]);
plot(t,y,'linewidth',1.5);
hold on;
plot(t,up,t,lo,'color','b');
plot(t,u + u_offset,'k','linewidth',1.5);
plot(t,up_u + u_offset,t,lo_u + u_offset,'color','k');
hold off;
grid on;

aa = subplot(2,2,2);
semilogx(log(omega),mag2db(up./up_u));
grid on;
axis([1e-2 2.5 -35 17]);

bb = subplot(2,2,4);
bodemag(G);
grid on;
% linkaxes([aa,bb],'xy');
axis([1e-2 21.5 -35 17]);
set(gca,'XTick',[]);



