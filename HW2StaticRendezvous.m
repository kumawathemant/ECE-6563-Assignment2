clear all, close all, clc

N=15; % Number of agents
dt=0.05; % numerical steplength
Tf=30; % final time
max_iter = Tf/dt;
delta=0.02; % termination condition

% Set up the networks (adjacency matrices)

% Complete Graph
AK=ones(N,N)-diag(ones(N,1));
LK=(N-1).*diag(ones(N,1))-AK;
eL=eig(LK); l2K=eL(2);

% Line Graph
AL=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
LL=2.*diag(ones(N,1))-AL; LL(1,1)=1; LL(N,N)=1;
eL=eig(LL); l2L=eL(2);

% Cycle Graph
AC=AL; AC(1,N)=1; AC(N,1)=1;
LC=2.*diag(ones(N,1))-AC;
eL=eig(LC); l2C=eL(2);


% Generate the initial conditions randomly
X=2.*rand(2,N);
% initial centroid
Xbar0 = 1/N*ones(1,N)*X';



% A is the adjacency matrix associated with the system
% Pick your favorite 
A=AL;
  
% Check for maximally allowable eps in Perron matrix
upper_bound=1/(max(sum(A)));

% Initialize plot
figure(1), hold on
set(gca,'xcolor','none','ycolor','none')
agenPlot = plot(X(1,:),X(2,:),'o','markersize',8,'MarkerFaceColor',[0.12,0.49,0.65],'MarkerEdgeColor','none');
axis([-0,2,-0,2])
  
for k = 1:max_iter-1
    

    u=zeros(2,N); %% Here is where we store the derivatives
    % FILL THIS PART!!!
    for i= 1:N 
        for j= 1:N 
            if (A(i,j) == 1)
                u(:,i) = u(:,i) - ( X(:,i) - X(:,j))*(norm(X(:,i) - X(:,j))-0.5);
                dist = norm(X(:,i) - X(:,j))

            end
        end
    end

    % Update the states using an Euler approximation
    for i=1:N
        X(:,i)=X(:,i)+dt.*u(:,i);
    end


    % Check if we have terminated
    xmax=0;
    for i=1:N
       xmax=max(xmax,norm(X(:,i)-Xbar0));
    end

    if xmax<delta
        break
    end

    % Update plot
    set(agenPlot,'xdata',X(1,:),'ydata',X(2,:))
    drawnow

end
