%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% RANSAC 1D - Multiple Target %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------- %
%          ALGORITHM OUTLINE          %
% ----------------------------------- %
% 0.1 Load measurements               %
% 0.2 Build constant matrices         %
% 1.  Make measurements into objects  %
%       (where each measurement has a %
%       .data,.t_step,.association)   %
% 2.  Build models for all measurement%
%       -s based on current time step %
%       measurements into a model     %
%       object vector                 %
% 3.0 Get rid of models that have the %
%       same initial and final inliers%
% 3.1 Choose model objects from vector%
%       that have a predefined number %
%       of associated measurements    %
% 4.  Apply KALMAN SMOOTHING to these %
%       models                        %
% 5. Delete associated measurements   %
%       from these models             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;
%% 0.1 Load measurements
data_file=load('1d_multTarg_ransac_data.mat');
t=data_file.t;      % Time array
dt=data_file.dt;    % Time step

%% 0.2 Build constant matrices
% Define Constants
N = 5; % measurement window length
n = 2;  % # of states (2 because of position and velocity)
m = 1;  % # of states I'm measuring (1 becauseI'm just obtaining position data)

% Dynamic model definition for 1D model
x_k=zeros(2,N);     % Initial States
trueModel = data_file.x;    % True state of target
y_k=data_file.y;    % Measurements: Columns 1-3 are true measurements. all other cols are from random noisy data
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

% Noise propagating portion;
noise_prop=inv(O'*inv(E)*O)*O'*inv(E); % size=nxN

%% 1.  Make measurements into objects %
%       (where each measurement has a %
%       .data,.t_step,.association)   %
meas = struct('data', cell(N, size(y_k,2)),...
    't_step', cell(N, size(y_k,2)),...
    'associated', cell(N, size(y_k,2)));
for kk=1:N
    for jj=1:size(y_k,2)
        meas(kk,jj).data = y_k(kk,jj);
        meas(kk,jj).t_step = kk;
        meas(kk,jj).associated = 0; % 0 = False
    end
end

%% 2. Build models for all measurement%
%       -s based on current time step %
%       measurements into a model     %
%       object vector                 %
tauR=3;
%Build model array (M) based on each current data points
M=struct('model',cell((N-1)*size(meas,2)*size(meas(end,:),2),1),...
    'assocmeas',cell((N-1)*size(meas,2)*size(meas(end,:),2),1));
modelnum=1;
for currdp=1:size(meas(end,:),2)
    y_currmeas=meas(end,currdp).data;
    % Build a model with currdp with respect to each data point at each
    % time step but the current one
    for ii=1:(size(meas,1)-1)
        for jj=1:size(meas,2)
            currInliers=0;
            tStep1=meas(ii,jj).t_step;
            y1=meas(ii,jj).data;
            x=noise_prop(:,[tStep1,meas(end,currdp).t_step])*[y1;y_currmeas];
            M(modelnum).model=[M(modelnum).model,x];
            for qq=1:(N-1) %time steps
                for rr=1:size(meas,2) %number of measurements per time step
                    diff=abs((x(1)+x(2)*qq)-meas(qq,rr).data);
                    if(diff<=tauR)
                        M(modelnum).assocmeas=[M(modelnum).assocmeas,meas(qq,rr)];
                        currInliers=currInliers+1;
                    end
                end
            end
            modelnum=modelnum+1;
        end    
    end
end

%% 3.0 Get rid of models that have the same initial and final inliers
eq_models=[];
currMod=1;
while(currMod<=size(M,1))
    if(size([M(currMod).assocmeas])>=1)
        [tstep_init, idx_init]=min([M(currMod).assocmeas.t_step]);
        inl_init_val=[M(currMod).assocmeas(idx_init).data];
        [tstep_fin, idx_fin] =max([M(currMod).assocmeas.t_step]);  
        inl_fin_val=[M(currMod).assocmeas(idx_fin).data];
        chckMod=1;
        while(chckMod<=size(M,1))
            %%Check inliers
            if(size([M(chckMod).assocmeas])>=1)
                [temp_tstep_init, temp_idx_init]=min([M(chckMod).assocmeas.t_step]);
                temp_inl_init_val=[M(chckMod).assocmeas(temp_idx_init).data];
                [temp_tstep_fin, temp_idx_fin] =max([M(chckMod).assocmeas.t_step]);  
                temp_inl_fin_val=[M(chckMod).assocmeas(temp_idx_fin).data];
                if((temp_tstep_init==tstep_init)&&(temp_inl_init_val==inl_init_val)&&...
                        (temp_tstep_fin==tstep_fin)&&(temp_inl_fin_val==inl_fin_val))
                    idx=find(eq_models==chckMod);
                    if(size(idx,2)==0)
                        M(chckMod)=[];
                        eq_models(end+1)=chckMod;
                    end
                end
            end
            chckMod=chckMod+1;
        end
    end
    currMod=currMod+1;
end

%% 3.1 Choose model objects from vector%
%       that have a predefined number %
%       of associated measurements    %
best_model.model=[];
best_model.assocmeas=[];
numbestmodels=0;
for ii=1:size(M,1)
    if (size(M(ii).assocmeas,2)>(N/6))
        numbestmodels=numbestmodels+1;
        best_model(numbestmodels).model=M(ii).model;
        best_model(numbestmodels).assocmeas=M(ii).assocmeas;
    end
end
%% Plot model on top of noisy data
model=zeros(size([best_model.model],2),size(t,2));
for ii=1:size([best_model.model],2)
    curr_best_model=best_model(ii).model;
    model(ii,:)=curr_best_model(1)+curr_best_model(2)*t;
end
clf;
figure(1);
for ii=1:size(model,1)
    plot(t(1,1:N),model(ii,1:N),'k--'); hold on;
end
plot(t(1,1:N),y_k(1:N,1),'b*'); hold on;
plot(t(1,1:N),y_k(1:N,2),'g*'); hold on;
%plot(t(1,1:N),y_k(1:N,3),'m.'); hold on;
%plot(t(1,1:N),y_k(1:N,4),'r.'); hold off;
%legend('Model','True data','Noise');
%pause(5)
%savefig('1d_1targ_ransac_comparison');
%close(gcf);
%% 4. Apply KALMAN SMOOTHING to these %
%       models                        %
P=inv(O'*inv(E)*O); % Error of covariance from paper
x_hat=zeros(2,N,size([best_model.model],2));
for ii=1:size([best_model.model],2) %Loop for each model
    curr_best_model=best_model(ii).model;
    x_hat(:,:,ii)=[curr_best_model(1)+curr_best_model(2)*t(1,1:N);curr_best_model(2)*ones(1,N)];
    for k=1:N %Loop through each time step (all time steps in time horizon) 
        if k~=1
            x_hat(:,k,ii)=A*x_hat(:,k-1,ii);
            P_new=A*P*A'+Qk;
        else
            P_new=P;
        end
        mod_inl=find([best_model(ii).assocmeas.t_step]==k);
        for j=1:size(mod_inl,2) %Each measurement inside model @ curr t_step
            L=P_new*C'*inv(R+C*P_new*C');
            P_new=(eye(n)-L*C)*P_new;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %x_hat(:,k,ii)=x_hat(:,k,ii)+L*(y_k(k,j)-C*x_hat(:,k,ii));
            x_hat(:,k,ii)=x_hat(:,k,ii)+L*(...
                best_model(ii).assocmeas(mod_inl(j)).data-C*x_hat(:,k,ii));
        end
    end
end

%% Plot smoothed model on top of noisy data
%clf;
figure(2);
for ii=1:size(x_hat,3)
    plot(t(1,1:N),x_hat(1,1:N,ii),'k--'); hold on;
end
plot(t(1,1:N),y_k(1:N,1),'b*'); hold on;
plot(t(1,1:N),y_k(1:N,2),'g*'); hold on;
%plot(t(1,1:N),y_k(1:N,3),'m.'); hold on;
%plot(t(1,1:N),y_k(1:N,4),'r.'); hold off;
%legend('Model','True data','Noise');
%pause(5)
savefig('1d_1targ_ransac_spread_smoothed');
%close(gcf);

%% 5. Delete associated measurements  %
%       from these models             %

