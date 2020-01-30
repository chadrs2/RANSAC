%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% RANSAC 1D - Multiple Target %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;
data_file=load('1d_multTarg_ransac_data.mat');
t=data_file.t;      % Time array
dt=data_file.dt;    % Time step

% Define Constants
N = 500; % measurement window length
n = 2;  % # of states (2 because of position and velocity)
m = 1;  % # of states I'm measuring (1 becauseI'm just obtaining position data)

% Dynamic model definition for 1D model
x_k=zeros(2,N);     % Initial States
trueModel = data_file.x;    % True state of target
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

% Solve initial states: x_k=inv(O'*inv(E)*O)*O'*inv(E)*y_k;
noise_prop=inv(O'*inv(E)*O)*O'*inv(E); % size=nxN
tauR=5;
numModels=0;
arrBestModels=[];
while(1)
    idxOfInliers=[];
    best_model=0;
    bestInliers=0;
    for l=1:1000
        currInliers=0;
        tStep1=randi(N);
        tStep2=randi(N);
        meas1=randi(size(y_k,2));
        meas2=randi(size(y_k,2));
        y1=y_k(tStep1,meas1);
        y2=y_k(tStep2,meas2);
        if(tStep1~=tStep2)
            currModel=noise_prop(:,[tStep1,tStep2])*[y1;y2];
            for ii=1:N
                for jj=1:size(y_k,2)
                    diff=abs((currModel(1)+currModel(2)*ii)-y_k(ii,jj));
                    if(diff<=tauR)
                        idxOfInliers=[idxOfInliers;[meas1, tStep1; meas2, tStep2]];
                        currInliers=currInliers+1;
                    end
                end
            end
            if(currInliers>=bestInliers)
                bestInliers=currInliers;
                best_model=currModel;
            else
                idxOfInliers=[];
            end
        end
    end
    best_model;
    bestInliers;
    if(bestInliers<(N))
        break;
    else
        for ii=1:size(idxOfInliers,1)*size(idxOfInliers,2)
            y_k(idxOfInliers(ii,2),idxOfInliers(ii,1))=0;
        end
        arrBestModels=[arrBestModels,best_model];
        numModels=numModels+1;
    end
end
% %% Plot model on top of noisy data
% model=best_model(1)+best_model(2)*t;
% clf;
% figure(1);
% plot(t(1,1:N),model(1,1:N),'k--'); hold on;
% plot(t(1,1:N),y_k(1:N,1),'b*'); hold on;
% plot(t(1,2:N),y_k(2:N,2),'r*'); hold off;
% legend('Model','True data','Noise');
% pause(5)
% savefig('1d_1targ_ransac_comparison');
% close(gcf);
% 
% %% Kalman Smoothing
% P=inv(O'*inv(E)*O); % Error of covariance from paper
% %x_hat=model;
% x_hat=[best_model(1)+best_model(2)*t(1,1:N);best_model(2)*ones(1,N)];
% for k=1:N
%     if k~=1
%         x_hat(:,k)=A*x_hat(:,k-1);
%         P_new=A*P*A'+Qk;
%     else
%         P_new=P;
%     end
%     for j=1:N*(m+1)
%         L=P_new*C'*inv(R+C*P_new*C');
%         P_new=(eye(n)-L*C)*P_new;
%         x_hat(:,k)=x_hat(:,k)+L*(y_k(k,randi(1))-C*x_hat(:,k));
%     end
% end
% 
% %% Plot smoothed model on top of noisy data
% clf;
% figure(1);
% plot(t(1,1:N),x_hat(1,1:N),'k--'); hold on;
% plot(t(1,1:N),y_k(1:N,1),'b*'); hold on;
% plot(t(1,2:N),y_k(2:N,2),'r*'); hold off;
% legend('Model','True data','Noise');
% pause(5)
% savefig('1d_1targ_ransac_smoothed');
% %close(gcf);