% 2019862s
% Yana Staneva

% Finite Volume Method 
% Specify value of m
m=2;

dx=0.01; % Define spatial step - if set to 0.001, remember to clear workspace
dt=dx; % Assign same temporal step, CFL satisfied

L = 10/m; % Space boundary parameter L
x = -L:dx:L; % Spatial domain

breakingTime = 2/m; % Breaking time
t = 0:dt:breakingTime; % Time interval

nu=dt/dx; % Courant number
n=length(x); % Number of spatial grid points
nt=length(t); % Number of temporal grid points

sechFunction = sech(m*x); % Hyperbolic secant function with m in the argument
u = zeros(1,n); % Initialize vector for iterations
uFiniteVolumeIter=sechFunction; % Initial condition for u
uFiniteVolume=zeros(nt, n); % Set 0 at the boundaries
uFiniteVolume(1,2:n-1)=sechFunction(2:n-1); % Initialize vector to store values of u

% Finite Volume Method
for k=1:nt
    
    fMinusSpace = sechFunction(1:n-2); % F_{j-1}^{n}
    fSpace = sechFunction(2:n-1); % F_{j}^{n}
    fPlusSpace = sechFunction(3:n); % F_{j+1}^{n}
    
    aHalfPlusSpace =0.5*(fPlusSpace.^2-fSpace.^2)./(uFiniteVolumeIter(3:n)-uFiniteVolumeIter(2:n-1)); % A_{j+1/2}^{n}
    aHalfMinusSpace =0.5*(fSpace.^2-fMinusSpace.^2)./(uFiniteVolumeIter(2:n-1)-uFiniteVolumeIter(1:n-2)); % A_{j-1/2}^{n}
   
    fHalfPlusSpace = 0.5.*((1+sign(aHalfPlusSpace)).*fSpace.^2 + (1-sign(aHalfPlusSpace)).*fPlusSpace.^2); % F_{j+1/2}^{n}
    fHalfMinusSpace = 0.5.*((1+sign(aHalfMinusSpace)).*fMinusSpace.^2 + (1-sign(aHalfMinusSpace)).*fSpace.^2); % F_{j-1/2}^{n}
    
    % Evaluate U_{j}^{n}
    u(2:n-1) = uFiniteVolumeIter(2:n-1)-0.5.*nu.*(fHalfPlusSpace - fHalfMinusSpace);
    
    % Store solution for next iteration
    uFiniteVolumeIter=u;
    % Update for the next iteration
    sechFunction=uFiniteVolumeIter;
    % Store output
    uFiniteVolume(k+1,:)=u;
end



finiteVolumeX=figure;
    plot(x(1:n),uFiniteVolume(26,:),'-k');    
    hold on
    plot(x(1:n),uFiniteVolume(51,:),':k');
    hold on
    plot(x(1:n),uFiniteVolume(76,:),'-.k');
    hold on
    plot(x(1:n),uFiniteVolume(101,:),'--k');
    title('Plot of u(x, t) as a function of spxace for \Delta x = 0.01 Using Finite Volume Method.');
    xlabel('x');
    ylabel('u(x,t)');
    legend('t=0.25','t=0.5','t=0.75','t=1','Location','northeast');
    %print(finiteVolumeX,'/Users/yanastaneva/Library/Mobile Documents/com~apple~CloudDocs/Year 5/Numerical Methods/Assignments/AC2/report/finiteVolumeXtest.jpeg','-djpeg');
    figure(finiteVolumeX);



