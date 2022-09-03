function [ R0,R1,C1,R2,C2,C0,gBest,Fitness_gBest ] = pso( V )
%% ----------------------------------------
% ----------- Scenario setup ------------
%----------------------------------------
load('Vmeasure.mat');
global Vbat0  
global PULSE
%PULSE =[0.84 1.39 2.1 4.2 8.4 12.6];(23&10'C)
%PULSE =[1.34 2.01 4.019 8.038 12.057];(45'C)
     Vbat0 = V(1:756012);
     PULSE =[0.84 1.39 2.1 4.2 8.4 12.6];
%% --------------------------------------
% ----------- PSO Parameters ------------
%----------------------------------------

NE              =  6;    % number of elements
N               =  100;  % number of particles

iter            =  100;   % number of iterations
iterCount       =  0;     % counter of iterations
c1              =  1;
c2              =  2;
w_max           =  1;
w_min           =  0.1;

%  -------------------------------------------
%  ------------ Boundary Settings ------------
%  -------------------------------------------
P_Max           =  [0.1 0.1 1000 0.1 5000 50000];     % [ max search range ]
P_Min           =  [10^-6 0 0 0 0 0];                 % [ min search range]
V_Max           =  (P_Max-P_Min+1)./80;         % [ max step ]
V_Min           =  - (P_Max-P_Min+1)./80;       % [ min step ]

P               =  zeros( N , NE );
V               =  zeros( N , NE );
pBest           =  zeros( N , NE );
gBest           =  zeros( iter , NE );
Fitness         =  zeros( N , 1 );
Fitness_pBest   =  zeros( N , 1 );
Fitness_gBest   =  zeros( iter , 1 );

%-----------------------------------------
% ------------ Initialization-------------
%-----------------------------------------
     for i = 1 : N               
         for j = 1 : NE
             P(i,j) =  (P_Min(1,j) + (P_Max(1,j)-P_Min(1,j))*rand());
             V(i,j) =  (V_Min(1,j) + (V_Max(1,j)-V_Min(1,j))*rand());
         end
     end
     
 while ( iterCount < iter )
%-------------------------------------
% ------------ Evaluation ------------
%-------------------------------------
 
   for i = 1 : N  
       Fitness(i,1) = Evaluation ( P(i,:) );   
   end

% ----------------------------------------------------------------------
%               Update pBest and pBest location & Fitness value
% ----------------------------------------------------------------------
    if iterCount == 0
        pBest = P;                                   % initial pBset
        Fitness_pBest = Fitness;
        [value, location] = min(Fitness_pBest);
        gBest( iterCount+1 , : ) = P(location,:);    % initial gBest
        Fitness_gBest( iterCount+1 , 1 ) = value;         
    else
         for j = 1 : N                               % update pBest
               if Fitness_pBest(j,1) >= Fitness(j,1) 
                  Fitness_pBest(j,1) = Fitness(j,1);
                  pBest(j,:) = P(j,:);
               end
          end
             [value, location] = min(Fitness_pBest);
             gBest(iterCount+1,:) = pBest(location,:); % update gBest
             Fitness_gBest(iterCount+1,1) = value;
    end
% ---------------------------------------------------------
%              Calculation velocities 
% ---------------------------------------------------------%

    w = w_max-(w_max-w_min)/iter*(iterCount+1);
    
    for i = 1 : N
        
        r1  =  rand( 1 , NE );
        r2  =  rand( 1 , NE );
        V(i,:) = w.*V(i,:) + r1.*c1.*(pBest(i,:)-P(i,:)) + r2.*c2.*(gBest(iterCount+1,:)-P(i,:));
        
                
        for j = 1 : NE
              if V(i,j) > V_Max(1,j)
                 V(i,j) = V_Max(1,j);
              end
              if V(i,j) < V_Min(1,j)
                 V(i,j) = V_Min(1,j);
              end
        end          

     end
 
% ---------------------------------------------------
%               Update locations
% ---------------------------------------------------
    for i = 1 : N

        P(i,:) = P(i,:) + V(i,:);
        P(i,:) =  P(i,:);
   
        for j = 1 : NE
              if P(i,j) > P_Max(1,j)
                 P(i,j) = P_Max(1,j);
              end
              if P(i,j) < P_Min(1,j)
                 P(i,j) = P_Min(1,j);
              end
        end
    end   
    
    

% ---------------------------------------------------
%               Next Iteration
% ---------------------------------------------------
disp(['iteration: ',num2str(iterCount+1)]);
disp(['best value ',num2str(Fitness_gBest(iterCount+1,1))]);
disp(['best setting',num2str( gBest(iterCount+1,:))]);

iterCount = iterCount + 1;    

    
 end
%%
figure;
plot(Fitness_gBest,'--.r');
hold on ;
xlabel('generation');
ylabel('fitness');
title('convergence history');
axis([0,100,1.5*10^(-3),10^(-2)]); 
grid

%%
R0 = gBest(100,1);
R1 = gBest(100,2);
C1 = gBest(100,3);
R2 = gBest(100,4);
C2 = gBest(100,5);
C0 = gBest(100,6);
end

     