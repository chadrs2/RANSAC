%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Code to Produce 1D Model of Data   %%
%   with multiple targets              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;

% Define constants
v = [randi(30);...
    %-randi(30);...
    randi(30)];  %m/s for all time
dt = .01;       %s
N = 1000;       % Number of time steps
M = 0;          % Number of false measuremnts per time step
t = (1:N)*dt;   % Time Array
Qk = 0.5;       % Process noise covariance
R = 0.1;        % Measurement noise covariance

% Position of targets as a function of time
x = zeros(N,2,size(v,1));   % 1st dim = num of time steps
                            % 2nd dim = num of states (pos, vel)
                            % 3rd dim = num of targets

% y -> contains the measurments, where the first
% three columns are the true measurement and the other 
% columns contain random noise measurements @ each
% time step
y = zeros(N,size(v,1)+M);

% Create position data
for targ=1:size(v,1)
    for k=2:N
        x(k,2,targ) = v(targ) + sqrt(Qk)*randn(1,1);
        x(k,1,targ)=x(k-1,1,targ)+dt*x(k,2,targ);
    end
end
x(:,1,1)=x(:,1,1)+1;

% Create noisy measurements
[~,idx]=max(abs(v));
for k=2:N
    for targ=1:size(v,1)
        % Add noise to true data
        y(k,targ) = x(k,1,targ)+sqrt(R)*randn(1,1);
    end
    for ii=1:M
        % Add random noise that are outliers
        %y(k,end)=max(max(max(abs(x(:,1,idx)))))*(-1+2.*rand(1,1));
    end
end
y(:,1)=y(:,1)+1;   
% Plot data
figure(1);clf;
plot(t,y(:,1),'b.')
hold on;
plot(t,y(:,2),'k.')
%hold on;
%plot(t,y(:,3),'g.')
%hold on;
%plot(t,y(:,end),'r.')
hold off;
xlabel('Time (s)');
ylabel('Position (m)');

% Save data
save('1d_multTarg_ransac_data.mat','t','v','x','y','dt','Qk','R')
savefig('1d_multTarg_ransac_data');
pause(2)
close(gcf);
