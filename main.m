%% IMPORTANT NOTES
% You must have YALMIP and OPTI toolbox to run this script.
% We've already prepared the conversion of PWA to the MLD model.
% The MLD model is saved in 'Threetank.m'
% But if you needed to change the model and convert the model again,
% follow these steps:

%% Conversion of PWA to MLD using HYSDEL
% HYSDEL 3.0 uses its own YALMIP toolbox, which is not compatible with the
% New version of YALMIP. So to make HYSDEL work, you should
% - Remove any version of the YALMIP toolbox from the MATLAB path.
% - Add HYSDEL 3.0 toolbox to MATLAB path.
% - Run the 2-line code below to extract the MLD model.

% clear classes
% hysdel3('Threetank.hys')

% - Remove the HYSDEL toolbox from the path. 
% - Add the new YALMIP toolbox to the MATLAB path.
% - Run the script 
% Feel free to ask your questions (ali.safi.eng@gmail.com)

clear;close all;clc;
yalmip('clear')
clear;clear classes;

%% Model data
Threetank
AffineModel
% MPC data
Q = 1*eye(S.nxr);
% Q(3,3)=100;
R = 1*eye(S.nur);
Rb = 0.001*eye(S.nub);
x_des=[.444;.35;.2];
x_0=[0;0;0];
N=2;
x=sdpvar(repmat(S.nxr,1,N+1),ones(1,N+1)); % solving the optimization problem for different initial conditions(consider it as sdpvar in the first element)
u=sdpvar(repmat(S.nur,1,N),ones(1,N));
ub=binvar(repmat(S.nub,1,N),ones(1,N));
z=sdpvar(repmat(S.nz,1,N),ones(1,N));
d=binvar(repmat(S.nd,1,N),ones(1,N));

constraints=[];
objective=norm(Q*(x{N}-x_des),1);

for i=1:N-1
    objective=objective+norm(Q*(x{i}-x_des),1) + norm(R*u{i},1) + norm(Rb*ub{i},1);
    
    constraints = [constraints, x{i+1} == S.A*x{i} + S.Bu(:,S.j.ur)*u{i} + S.Bu(:,S.j.ub)*ub{i} + S.Baux(:,S.j.z)*z{i} + S.Baux(:,S.j.d)*d{i} + S.Baff,...
        S.Ex*x{i} + S.Eu(:,S.j.ur)*u{i} + S.Eu(:,S.j.ub)*ub{i} + S.Eaux(:,S.j.z)*z{i} + S.Eaux(:,S.j.d)*d{i}<=S.Eaff];
    
end

controller = optimizer(constraints, objective , sdpsettings('Solver','cbc','verbose',0),x{1},[u{1};ub{1};x{2};z{1};d{1}]);

%% closed-loop simulation
horizon=60;
x=zeros(S.nxr,horizon+1);
u=zeros(S.nur,horizon);
ub=zeros(S.nub,horizon);
z=zeros(S.nz,horizon);
d=zeros(S.nd,horizon);
x(:,1)=x_0;

for i=1:horizon
    for q=1:size(x,1)
        if(x(q,i)<0)
            x(q,i)=0;
        end
    end
    f=controller(x(:,i)+0.03*rand(1));
    u(:,i)=f(1:S.nur);
    ub(:,i)=f(S.nur+1:S.nur+S.nub);
    if ub(1,i) && ub(2,i)
        j=1;
    elseif ub(1,i) && ~ub(2,i)
        j=2;
    elseif ~ub(1,i) && ub(2,i)
        j=3;
    elseif ~ub(1,i) && ~ub(2,i)
        j=4;
    end
    
x(:,i+1)=S.model(j).A*x(:,i)+S.model(j).B*u(:,i)+S.model(j).d;
end

%%
mpcawn.x=x;
mpcawn.u=u;
mpcawn.ub=ub;
figure('name','states')
subplot(3,1,1)
plot(x(1,:),'linewidth',1.5); axis tight
xlabel('sample')
ylabel('water level 1')
hold on
plot([0 horizon],[x_des(1) x_des(1)],'linewidth',1.5)
legend(["x_1","desired"])


subplot(3,1,2)
plot(x(2,:),'linewidth',1.5); axis tight
xlabel('sample')
ylabel('water level 2')
hold on
plot([0 horizon],[x_des(2) x_des(2)],'linewidth',1.5)
legend(["x_2","desired"])


subplot(3,1,3)
plot(x(3,:),'linewidth',1.5); axis tight
xlabel('sample')
ylabel('water level 3')
hold on
plot([0 horizon],[x_des(3) x_des(3)],'linewidth',1.5)
legend(["x_3","desired"])



figure('name','flow')
stairs(u','linewidth',1.5); axis tight
legend(["valve 1";"valve 2"])
grid on

figure('name','switching')
stairs(ub','linewidth',1.5); axis tight
legend(["valve 1";"valve 2"])
grid on