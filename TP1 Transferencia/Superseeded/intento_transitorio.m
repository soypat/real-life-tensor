%% 1-D FVM code. Bar with conduction
Nx=40; % Number of volumes
L=1; %metro
x = linspace(0,L,Nx+1); % Surface positions for volumes
dx = L/Nx;
xmid = 0.5*(x(1:Nx) + x(2:Nx+1)); %Volume midpoints

%Let's solve a transient heat conduction 1-D problem.
%simple!
%dT/dt - alpha*div(T) = 0
%For our problem. 1D: dT/dt=alpha*(d^2T/dx^2)
%alpha=k/(rho*cp)

% Initial conditions: Temperature along bar (T=25 constant)
% Border conditions: T(1)=0
%                    T(Nx)=30
%                    
%% Bar dimensions
A=0.1; %m^2 Area
% L=L; %Meters, declared above
cp=510; % Joule per kilogram kelvin (steel)
k=43;
rho=9750; % kilogram per meter cubed
alpha=k/rho/cp;
%% Solver
T=25*ones(Nx,1);
T0=0;%Border condition
TN=30;
Tbc = [T0; T ; TN]; %We have two border conditions
dt=200; %seconds
tfinal=dt*2000;
t=0;%Initialize start time 
%{
T0     T(1)                T(Nx)     TN
|-------|------- ...  -------|-------|
     <dT[1]>             <dT[Nx]>


BAR VOLUME
               
            |-----------|
            |           |
EAST (Cond.)|           | WEST (Conduction)
            |           |
            |           |
            |-----------|
%}
vol=dx*A;
Q=zeros(Nx,1);
while (t<tfinal)
    Tbc = [T0; T ; TN];
    for i = 1:Nx
        dTdxwest = Tbc(i+1)-Tbc(i);
        dTdxeast = Tbc(i+1)-Tbc(i+2);
        dTdt= -alpha*A/vol*(dTdxwest+dTdxeast);
        Q(i)=dTdt*dt;
    end
    t=t+dt;
    T=T+Q;
    %Plot current solution
     stairs(x,[T; T(Nx)]);
     axis([0, 1, min(T0,TN), max(T0,TN)]);
     grid on;
     drawnow;
end
