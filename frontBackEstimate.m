% 2019862s
% Yana Staneva

% Estimate Back and Front of Wave

dx=0.01; % Define spatial step - if set to 0.001, remember to clear workspace
dt=dx; % Assign same temporal step, CFL satisfied

X=2; % Spatial boundary
x=0:dx:X; % Spatial domain
xb=0.1; % Initial condition for back of the wave
xf=0.3; % Initial condition for front of the wave
Tf=2; % Temporal boundary
t = (0:dt:2); % Temporal domain

n=floor(X/dx); % Number of grid points in space
n1=floor(xb/dx); % Grid points before the back of the wave at t=0
n2=floor(xf/dx); % Grid points before the front of the wave at t=0
uExactIter = zeros(1,n+1); % Initialize exact solution
uExact = zeros(1,n+1); % Initialize variable to store the exact solution

nt=floor(Tf/dt); % Number of grid points in time domain
nu=dt/dx; % Courant number

for i = 1:nt+1
    polyBack = [1 -0.1 1 -t(i)-0.1]; % Polynomial interpolant for back
    polyFront = [1 -0.3 1 -t(i)-0.3]; % Polynomial interpolant for front
  
    xback = roots(polyBack); % Compute the roots of the back polynomial
    xback = xback(imag(xback)==0); % Store only real roots
    xb(i)=xback; % Store the computed root for the current iteration
    
    xfront = roots(polyFront); % Compute the roots of the front polynomial
    xfront = xfront(imag(xfront)==0); % Store only real roots
    xf(i)=xfront; % Store the computed root for the current iteration
    
    for j = 1:n+1
         if x(j)>=xb(i) && x(j)<=xf(i) % For x between calculated xb and xf
            uExactIter(j) = 1; % Set the amplitude of the wave to be 1
         else 
            uExactIter(j) = 0; % Otherwise it is 0
         end
    end
    uExact(i,:) = uExactIter; % Save the obtained value to use in the next iteration
end



u0=[zeros(1,n1),ones(1,n2-n1+1),zeros(1,n-n2)]; % Initial condition for u
u = u0; % Initialize u to IC values
uNumerical(1,:)=u0; % Initialize vector to store the numerical solution
time(1,1)=0; % Initialize vector to store the temporal interval
uNumericalIter=u0; % Initialize vector to store the values for the next iteration
for k=1:nt
    t=k*dt; % Current time for the iteration

    % Velocity functions 
    a=(1+x.^2)./((1+x.^2).^2 + 2*x*t);
    aHalfPlusTime=(1+x.^2)./((1+x.^2).^2+2*x*t+x*dt); % A^{n+1}
    aHalfPlusSpace=(1+x.^2+x.*dx+0.25*dx.^2)./(((1+x.^2+x.*dx+0.25*dx.^2)).^2 + 2*x*t+dx*t); % A_{j+1/2}
    aHalfMinusSpace=(1+x.^2-x.*dx+0.25*dx.^2)./(((1+x.^2-x.*dx+0.25*dx.^2)).^2 + 2*x*t-dx*t); % A_{j-1/2}
    
    % Preallocate values of neighbouring U nodes
    uPlusSpace=(uNumericalIter(3:n+1)-uNumericalIter(2:n)); % U_{j+1}
    uMinusSpace=(uNumericalIter(2:n)-uNumericalIter(1:n-1)); % U_{j-1}
  
    % Evaluate U_{j}^{n}
    u(2:n) = uNumericalIter(2:n)-(0.5*aHalfPlusTime(2:n).*nu.*(uNumericalIter(3:n+1)-uNumericalIter(1:n-1)))+ ...
          (0.5*nu.^2.*a(2:n)).*(aHalfPlusSpace(2:n).*uPlusSpace-aHalfMinusSpace(2:n).* uMinusSpace);
    
    % Update time with temporal step
    time(k+1,1)=k*dt;
    % Store value of U for next iteration
    uNumericalIter=u;
    % Store Output
    uNumerical(k+1,:)=u;
end

% Specify error - a change in the error will not be automatically reflected
% in the name of the figure produced or the path stored, so these should be
% changed as well, see below
error = 0.5;

for l = 1:n,
    for m = 1:n,
        if abs(uNumerical(l,m))>error, % If |u|>error
           frontEstimate(l,:) = [time(l), x(m)]; % Estimate front
        end
    end
    for h = n:-1:1,
        if abs(uNumerical(l,h))>error, % If |u|>error
           backEstimate(l,:) = [time(l), x(h)]; % Estimate back
        end
    end
end


frontBackPrediction=figure;
    plot(frontEstimate(:,1),frontEstimate(:,2),'-k');
    hold on
    plot(backEstimate(:,1), backEstimate(:,2),':k');
    hold on
    plot(time, xf(:),'-.k');
    hold on
    plot(time, xb(:), '--k');
    title('Approximation of Back and Front of Wave Using Numerical Solution for \Delta x = 0.01 and \epsilon = 0.5.'); % Name should be changed 
    xlabel('t');
    ylabel('Back and Front of Wave.');
    legend('Front Numerical','Back Numerical','Front Exact','Back Exact','Location','northwest');
    %print(frontBackPrediction,'/Users/yanastaneva/Library/Mobile Documents/com~apple~CloudDocs/Year 5/Numerical Methods/Assignments/AC2/report/frontBackPrediction0001eps05test.jpeg','-djpeg');
    figure(frontBackPrediction)
    
