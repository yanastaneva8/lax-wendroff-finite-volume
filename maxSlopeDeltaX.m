% 2019862s
% Yana Staneva

% Calculate Maximum Slope for Different Spatial Steps


deltaX = 0.001:0.001:0.02; % Vector of spatial steps

for i=1:length(deltaX)
    m = 2; % Keep m = 2
    L = 10/m; % Space boundary parameter L
    
    x = -L:deltaX(i):L; % Spatial domain to be updated at every iteration
    n=length(x); % Number of spatial grid points
    
    dt=deltaX(i); % Extract i-th value of deltaX vector
    breakingTime = 2/m; % Breaking time
    t = 0:dt:breakingTime; % Time interval
    nt=length(t); % Number of temporal grid points
    
    nu=dt/deltaX(i); % Courant number updated as well 
    
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
        slope(j, :) = abs((uJplus1 - uJ))/deltaX(i); % Compute absolute value of difference 
    end
    
    slopeMatrix = slope'; % Transpose the slope matrix as max only works for columns
    maximumSlope = max(slopeMatrix); % Get largest value of slope
    maximumSlopeIter = maximumSlope(k); % Store the max value from iteration
    maxSlopeAtBreakTime(i) = maximumSlopeIter; % Update for next iteration 
    fprintf('%8.3f & %8.3f \\\\ \n', deltaX(i), maxSlopeAtBreakTime(i)) % Prints results in LaTeX table format
end
