%% Figure 1: Samantha's Exploration 

% define probability of transition
p11 = 0.3; p12 = 0.7; p22 = 0.5; p23 = 0.5; p31 = 0.2; p33 = 0.80;

% define stochastic matrix
P = [p11 p12 0; 0 p22 p23; p31 0 p33];

% prepare simulation
time = 30; % time step
location = zeros(1,time); % store Samantha's location

% start simulation
location(1) = 1; % Samantha starts at location 1

for t=2:time
    possibleNextLoc = P(location(t-1),:);
    nextLoc = sample_state(possibleNextLoc);
    location(t) = nextLoc;
end

plot(1:time, location, '.--','Markersize',20,'Linewidth',1.5)
hold on 
yline(1)
hold on 
yline(2)
hold on
yline(3)
yticks([1 2 3])
ylim([0.9 3.1])
xlabel('time (days)')
ylabel("Samantha's Location")

%% Figure 2: Samantha's stationary distribution

% define probability of transition
p11 = 0.3; p12 = 0.7; p22 = 0.5; p23 = 0.5; p31 = 0.2; p33 = 0.80;

% define stochastic matrix
P = [p11 p12 0; 0 p22 p23; p31 0 p33];

% define starting distribution 
mu0 = [0.4 0.5 0.1];

% prepare computation
time = 20; % time step
location = zeros(3,time); % store Samantha's location

% start computation
location(:,1) = mu0';

for t=2:time 
    distribution = location(:,t-1)'*P;
    location(:,t) = distribution';
end

% get theoretical result 
[V,D] = eig(P');
pi = V(:,find(diag(D) >= 1-1e-6));
pi = pi/sum(pi); % matlab normalized w.r.t 2-norm

plot(1:time, location,'Linewidth',1.5)
hold on
yline(pi(1),'--k','Linewidth',1)
hold
yline(pi(2),'--k','Linewidth',1)
hold
yline(pi(3),'--k','Linewidth',1)
hold
legend('Prob. at Field 1','Prob. at Field 2','Prob. at Field 3') 
xlabel('time (days)')
ylabel('probability')
ylim([0 1])
xlim([1 time])


%% Buffon's Needle Problem
repeat = 3;
trial = 300;
L = 1; % needle length

result = zeros(repeat,trial);

for r = 1:repeat
    count = 0; % just for computational complexity
    for t=1:trial
        X = rand * L/2; % generate the midpoint
        theta = rand * pi/2; % generate the angle of the needle
        intersect = (X < L*cos(theta)/2);

        if intersect
            count = count + 1;
        end

        result(r,t) = count/t;
    end
end


plot(1:trial, result, 'Linewidth',1.5)
hold on
yline(2/pi,'--k','Linewidth',1.5)
ylim([0 1])
xlabel('number of trials')
ylabel('probability of intersecting')

%% Tends towards low energy

time = 30000;
temp = [1.3, 2.7, 4];

for tt = 1:length(temp)
    beta = 1/temp(tt);

    n = 50; 
    N = n^2; % system size
    sigmaM = zeros(time,N); % Metropolis-Hastings
    sigmaG = zeros(time,N); % Glauber Dynamics
    energyM = zeros(time,1);
    energyG = zeros(time,1); 

    % generate lattice with closed boundary condition
    A = delsq(numgrid('S',n+2));
    G = graph(A,'omitselfloops');

    % add boundary conditions
    for i=1:n
        G = addedge(G,(i-1)*n+1,i*n,1);
        G = addedge(G,i,n*(n-1)+i,1);
    end
    H = adjacency(G); % the Hamiltonian is the adjacency matrix

    sigma0 = ((rand(1,N) > 0.5)-1/2)*2; % start randomly
    sigmaM(1,:) = sigma0;
    sigmaG(1,:) = sigma0;
    energyM(1) = -sigma0*H*sigma0';
    energyG(1) = -sigma0*H*sigma0';

    % start simulation
    for t=2:time
        % Metropolis-Hastings
        current = sigmaM(t-1,:);

        % generate new by flipping one site
        new = current;
        randSite = randi(N);
        new(randSite) = new(randSite)*-1;

        Hcurrent = -current*H*current'/2;
        Hnew = -new*H*new'/2;

        r = min(exp(-beta*(Hnew - Hcurrent)),1);

        if rand < r
            sigmaM(t,:) = new;
            energyM(t) = Hnew;
        else
            sigmaM(t,:) = current;
            energyM(t) = Hcurrent;
        end

        % Glauber Dynamics
        current = sigmaG(t-1,:);
        randSite = randi(N);

        S = sum(H(randSite,:).*current);
        probUp = (1 + tanh(beta*S))/2;

        if rand < probUp
            current(randSite) = 1;
        else
            current(randSite) = -1;
        end
        sigmaG(t,:) = current;
        energyG(t) = -current*H*current'/2;
    end

    subplot(2,length(temp),tt)
    plot(1:time,energyM,'Linewidth',1.5)
    hold on
    plot(1:time,energyG,'--','Linewidth',1.5)
    xlabel('iterations')
    ylabel('energy')
    title(['tempature =',num2str(temp(tt))])
    legend('Metropolis-Hastings','Glauber Dynamics')
    
    subplot(2,length(temp),tt+length(temp))
    sigmaMatrix = reshape(sigmaG(time,:),n,n);
    imagesc(sigmaMatrix)
    colormap(gray)
end



%% Function definitions

% draw samples for pmf
function x = sample_state(p)
cdf = cumsum(p);
xi = rand;
for i=1:length(p)
    if cdf(i) > xi
        x = i;
        break
    end
end
end


