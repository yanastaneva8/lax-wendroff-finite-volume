% 2019862s
% Yana Staneva

% Calculate the Maximum Slope of U as a function of Time for m=2,3,4,5.

mValues = [2, 3, 4, 5]; % Possible values for m
mGradients = zeros(length(mValues), 101); % Initialize vector for the slopes


for i=1:length(mValues)
   
m = mValues(i); % Take value from vector of possible m values

dx=0.01; % Define spatial step - if set to 0.001, remember to clear workspace
dt=dx; % Assign same temporal step, CFL satisfied

L = 10/m; % Space boundary parameter L
x = -L:dx:L; % Spatial domain
n=length(x); % Number of spatial grid points

breakingTime = 2/m; % Breaking time
t = 0:dt:breakingTime; % Time interval
nt=length(t); % Number of temporal grid points

nu=dt/dx; % Courant number

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

    timePoints=nt; 
    slope = zeros(nt, n-1); % Initialize matrix to store values of slope
    for j=1:timePoints % Iterate for all t
        uJ = uFiniteVolume(j,1:n-1); % Get U_{j}^{n}
        uJplus1 = uFiniteVolume(j, 2:n); % Get U_{j+1}^{n}
        slope(j, :) = abs((uJplus1 - uJ))/dx; % Compute absolute value of difference 
    end

    slopeMatrix = slope'; % Transpose the slope matrix as max only works for columns
    maximumSlope = max(slopeMatrix); % Get largest value of slope
    mGradients(i, 1:k) = maximumSlope; % Store the max value in a matrix

end    


maxSlope=figure;
    plot(0:dt:1, mGradients(1,:), '-k', ...
     0:dt:2/3, mGradients(2,1:length(0:dt:2/3)),':k', ...
     0:dt:1/2, mGradients(3,1:length(0:dt:1/2)), '-.k', ...
     0:dt:2/5, mGradients(4,1:length(0:dt:2/5)), '--k');
    title('Plot of the Maximal Slope of u(x,t) as a function of time for \Delta x = 0.01 in [0,t_{m}]');
    xlabel('t');
    ylabel('du/dx');
    legend('m=2','m=3','m=4','m=5','Location','northeast');
    %print(maxSlope,'/Users/yanastaneva/Library/Mobile Documents/com~apple~CloudDocs/Year 5/Numerical Methods/Assignments/AC2/report/maxSlopetest.jpeg','-djpeg');
    figure(maxSlope);

