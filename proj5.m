%% Project 5

load('proj5.mat')       % data loaded from text files

x01 = 0.1:0.1:0.9;          % dx = 1/10
x001 = 0.01:0.01:1-0.01;    % dx = 1/100

dt1 = 0.05;
dt2 = 5e-5;
t1 = 20;    % steps
t2 = 120;

u01 = analytical(x01,dt1,t1);       % dx = 0.1, dt = 0.05 at t1
u012 = analytical(x01,dt1,t2);      % dx = 0.1, dt = 0.05 at t2
u001 = analytical(x001,dt2,t1);     % dx = 0.01, dt = 5e-5 at t1
u0012 = analytical(x001,dt2,t2);    % dx = 0.01, dt = 5e-5 at t2

figure
plot(x01,u01)
figure
plot(x01,u012)
figure
plot(x001,u001)
figure
plot(x001,u0012)


%% one dim x = 1/10

figure  % at t1
hold on
plot(x01,onedimforward01(t1/5,:),'Displayname','FW Euler')
plot(x01,onedimbackward01(t1/5,:),'Displayname','BW Euler')
plot(x01,onedimcrank01(t1/5,:),'Displayname','Cranck-Nicolson')
plot(x01,u01./max(u01),'Displayname','Exact') % normalized
ylabel('Temperature diffusion','Fontsize',16)
xlabel('Length','Fontsize',16)
lgd = legend('show','Location','Northwest');
lgd.FontSize = 14;
%%
figure  % at t2
hold on
plot(x01,onedimforward01(t2/5,:),'Displayname','FW Euler')
plot(x01,onedimbackward01(t2/5,:),'Displayname','BW Euler')
plot(x01,onedimcrank01(t2/5,:),'Displayname','Cranck-Nicolson')
plot(x01,u012./max(u012),'Displayname','Exact')
ylabel('Temperature diffusion','Fontsize',16)
xlabel('Length','Fontsize',16)
lgd = legend('show','Location','Northwest');
lgd.FontSize = 14;

%% one dim x = 1/100

figure  % t1
hold on
plot(x001,onedimbackward001(4,1:99),'Displayname','BW Euler')
% plot(x001,onedimcrank001(20,:),'Displayname','Crank-Nicolson') %not right
plot(x001,u001./max(u001),'Displayname','Exact')
ylabel('Temperature diffusion','Fontsize',16)
xlabel('Length','Fontsize',16)
lgd = legend('show','Location','Northwest');
lgd.FontSize = 14;

%%
figure  % t2
hold on
plot(x001,onedimbackward001(24,1:99),'Displayname','BW Euler')
% plot(x001,onedimcrank001(24,:))
plot(x001,u001./max(u001),'Displayname','Exact')
ylabel('Temperature diffusion','Fontsize',16)
xlabel('Length','Fontsize',16)
lgd = legend('show','Location','Northwest');
lgd.FontSize = 14;


%% two dim

% % 10x10
% figure
% imagesc(twodim10)
% colormap hot
% axis xy
% xlabel('X','Fontsize',16)
% ylabel('Y','Fontsize',16)
% title('Temperature','Fontsize',14)

% dt = .0000625, bc explicit - stability condition: dt/dx^2 le 1/2
% so dt is determined from dx in the code

% 20x20
figure
imagesc(twodim20)
colormap hot
axis xy
xlabel('X','Fontsize',16)
ylabel('Y','Fontsize',16)
title('Temperature','Fontsize',14)

% heat from two sides
figure
imagesc(twodim202)
colormap hot
axis xy
xlabel('X','Fontsize',16)
ylabel('Y','Fontsize',16)
title('Temperature','Fontsize',14)

% 40x40
figure
imagesc(twodim40)
colormap hot
axis xy
xlabel('X','Fontsize',16)
ylabel('Y','Fontsize',16)
title('Temperature','Fontsize',14)


