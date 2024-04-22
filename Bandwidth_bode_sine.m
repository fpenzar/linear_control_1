clear all;
clc;

s = tf('s');
G = 25/(s^2 + 2*0.254*5*s + 25);

w = logspace(-1,1.2,1000);          % Define frequency range for Bode plot
[M,P,w_out] = bode(G,w);            % Get magnitude and phase of G
M = mag2db(squeeze(M));
P = squeeze(P);

i_BW = find(M <= M(1) - 3,1,'first');   % Get the index of the bandwidth frequency
omega_BW = w_out(i_BW);                 % Get omega_BW

% Bode plots
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultLineLineWidth', 2);

figure(1);
a1 = subplot(2,1,1);
semilogx(w_out,M,'Color',[84 130 15]./255);
hold on;
semilogx(w_out,(M(1) - 3)*ones(size(w_out)),'--b');
scatter(w_out(i_BW),M(i_BW),25,'markerFaceColor','k','MarkerEdgeColor',' k');
hold off;
grid on;
ylabel('$\vert G_{cl} \vert$ in dB');
set(gca,'xtick',[]);
legend({'$\vert G_{cl} \vert$' ,'$-3$dB'},'interpreter','latex');

b1 = subplot(2,1,2);
semilogx(w_out,P,'Color',[84 130 15]./255);
hold on;
semilogx(w_out,-180*ones(size(w_out)),'--k');
hold off;
grid on;
xlabel('$\omega$ in rad/s');
ylabel('$\angle G_{cl}$ in deg');
linkaxes([a1,b1],'x');
xlim([w_out(1) w_out(end)]);

%% Sine wave responses
t = linspace(0, 10, 10000);     % Time Vector
omega = [0.5 1 2.5 5 7.4 14];
for i = 1:length(omega)    
    u = sin(omega(i)*t);        % Input
    y = lsim(G, u, t);          % Calculate System Response
    figure(2);
    subplot(3,2,i)
    plot(t,u,'k');
    hold on;
    plot(t,y,'Color',[4 145 245]./255);
    hold off;
    grid on;
    if (i >= 5)
        xlabel('$\omega$ in rad/s');
    end
    if (mod(i,2) ~= 0)
        ylabel('$y(t)$');
    end
    ylim([-2 2]);
    title(strcat(['$r(t) = \sin(',num2str(omega(i)),'t)$']))
end


