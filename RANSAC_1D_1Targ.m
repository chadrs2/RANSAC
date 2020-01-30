%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% RANSAC 1D - 1 Target %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;
data_file=load('1d_ransac_data.mat');
t=data_file.t;      % Time array
dt=data_file.dt;    % Time step

% Define Constants
N = 500; % measurement window length
n = 2;  % # of states (2 because of position and velocity)
m = 1;  % # of states I'm measuring (1 becauseI'm just obtaining position data)

% Dynamic model definition for 1D model
x_k=zeros(2,N);     % Initial States
x = data_file.x;    % True state of target
y_k=data_file.y;    % Measurements: 1st column are measurements from the system. all ofther are from random noisy data
a = [0 1; 0 0];     % Continuous time system matrix
A = [1 dt; 0 1];    % Discretized system matrix
C = [1 0];          % Measurement matrix
Gk = (2*eye(2)+a*dt+a^2*dt^2/3+a^3*dt^3/12)*dt/2*[0;1]; % Measurements process noise matrix
Q = data_file.Qk;   % Continuous time process noise covariance
R = data_file.R;    % Continuous time measurement noise covariance
Qk= Gk*Q*Gk';       % State/Processing noise discretized

% Other matrix definitions
O=zeros(m*N,n);     % Matrix describing the expected measurement's evolution form initial state to subsequent time steps
G=zeros(m*N,n*N);   % Propagates the process noise 

% Build/Fill in Matrices
%---------O---------%
for i=1:N
    if i~=1
        O(i,:)=C*(A^(i-1));
    else
        O(i,:)=C;
    end
end
%---------G---------%
for i=1:m*N
    for j=1:n:n*N
        if (i==1 || j==1)
            G(i,j:j+n/2)=zeros(n,1);
        elseif (i==ceil(j/n))
            G(i,j:j+n/2)=C;
        elseif (i>j)
            G(i,j:j+n/2)=C*(A^(i-j));
        else
            G(i,j:j+n/2)=0;
        end
    end
end
%---------E---------%
Qblock=kron(eye(m*N),Qk);   % size=nNxnN
Rblock=kron(eye(m*N),R);    % size=mNxmN
E=G*Qblock*G'+Rblock;       % size=mNxmN

%% Solve initial states: x_k=inv(O'*inv(E)*O)*O'*inv(E)*y_k;
noise_prop=inv(O'*inv(E)*O)*O'*inv(E); % size=nxN
tauR=5;
best_model=0;
bestInliers=0;

for l=1:10000
    currInliers=0;
    tStep1=randi(10);
    tStep2=randi(10);
    y1=y_k(tStep1,randi(2));
    y2=y_k(tStep2,randi(2));
    if(tStep1~=tStep2)
        x=noise_prop(:,[tStep1,tStep2])*[y1;y2];
        for ii=1:N
            for jj=1:(m+1)
                diff=abs((x(1)+x(2)*ii)-y_k(ii,jj));
                if(diff<=tauR)
                    currInliers=currInliers+1;
                end
            end
        end
        if(currInliers>=bestInliers)
            bestInliers=currInliers;
            best_model=x;
        end
    end
end
best_model
bestInliers

%% Plot model on top of noisy data
model=best_model(1)+best_model(2)*t;
clf;
figure(1);
plot(t(1,1:N),model(1,1:N),'k--'); hold on;
plot(t(1,1:N),y_k(1:N,1),'b*'); hold on;
plot(t(1,2:N),y_k(2:N,2),'r*'); hold off;
legend('Model','True data','Noise');
pause(5)
savefig('1d_1targ_ransac_comparison');
close(gcf);

%% Kalman Smoothing
P=inv(O'*inv(E)*O); % Error of covariance from paper
%x_hat=model;
x_hat=[best_model(1)+best_model(2)*t(1,1:N);best_model(2)*ones(1,N)];
for k=1:N
    if k~=1
        x_hat(:,k)=A*x_hat(:,k-1);
        P_new=A*P*A'+Qk;
    else
        P_new=P;
    end
    for j=1:N*(m+1)
        L=P_new*C'*inv(R+C*P_new*C');
        P_new=(eye(n)-L*C)*P_new;
        x_hat(:,k)=x_hat(:,k)+L*(y_k(k,randi(1))-C*x_hat(:,k));
    end
end

%% Plot smoothed model on top of noisy data
clf;
figure(1);
plot(t(1,1:N),x_hat(1,1:N),'k--'); hold on;
plot(t(1,1:N),y_k(1:N,1),'b*'); hold on;
plot(t(1,2:N),y_k(2:N,2),'r*'); hold off;
legend('Model','True data','Noise');
pause(5)
savefig('1d_1targ_ransac_smoothed');
%close(gcf);