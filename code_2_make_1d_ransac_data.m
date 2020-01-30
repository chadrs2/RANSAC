%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Code to Produce 1D Model of Data   %%
% -------------------------------------%%
%       P(k)=P(k-1)+v(k)*dt+w(k)        %
%       ------------------------        %
% VARIABLE DEFINITIONS:                 %
% P - Position                          %
% v - velocity (constant)               %
% dt - time step (constant)             %
% N - number of data points/time steps) %
% w - noise, where wk~N(0,Qk)           %
%   > Qk - variance in noise            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;

% Define constants
v = 5;          %m/s for all time
dt = .01;       %s
N = 1000;       % Number of time steps
M = 1;          % Number of false measuremnts per time step
t = (1:N)*dt;   % Time Array
Qk = 0.5;       % Process noise covariance
R = 0.1;        % Measurement noise covariance

% Position of target as a function of time
x = zeros(N,2);

% y -> contains the measurments, where the first
% column is the true measurement and the other 
% columns contain random noise measurements @ each
% time step
y = zeros(N,1+M);

% Create position data
for k=2:N
    x(k,2) = v + sqrt(Qk)*randn(1,1);
    x(k,1)=x(k-1,1)+dt*x(k,2);
end

% Create noisy measurements
for k=2:N
    % Add noise to true data
    y(k,1) = x(k,1)+sqrt(R)*randn(1,1);
    
    for ii=1:M
        % Add random noise that are outliers
        y(k,2)=max(max(max(x(:,1))))*(-1+2.*rand(1,1));
    end
end
   
% Plot data
figure(1);clf;
plot(t,y(:,1),'b.')
hold on;
for ii = 1:M
    plot(t,y(:,1+ii),'r.')
end
xlabel('Time (s)');
ylabel('Position (m)');

% Save data
save('1d_ransac_data.mat','t','x','y','dt','Qk','R')
savefig('1d_ransac_data_rand');
%close(gcf);

