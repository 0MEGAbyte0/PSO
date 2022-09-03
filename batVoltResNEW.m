function [ maxError,Vbat0, Vbat  ] = batVoltResNEW( R0,R1,C1,R2,C2,C0,VSTART )
%% 2RC 參數設定
% 用來驗證擬合結果
% % 畫出給定參數的電壓響應

Qmax = 4.2*3600;
Umax = 4.2;
Umin = 2.5;

load('Vmeasure.mat');
    Vbat0 = V(1:756012);
    PulseStrenth = [0.84 1.39 2.1 4.2 8.4 12.6];
%PulseStrenth = [1.34 2.01 4.019 8.038 12.057]

%% EHPPC 參數設定
pusleStart = 0;    % [second]
pulseDura = 30.01;     % [second]
relaDura = 600;     % [second]
sampleTime = 0.01;  % [second]
%PulseStrenth =  12.6 ;
%PULSE =[0.84 1.39 2.1 4.2 8.4 12.6];

%% 初始化

StartID = pusleStart/sampleTime;    % data point index

% index by time [second]
PP = 0:sampleTime:pulseDura;        % Duration of pulse
NP = pulseDura/sampleTime;          % Number of Pulse Point
RP = 0:sampleTime:relaDura;         % Duration of relaxation
NR = relaDura/sampleTime;           % Number of Relaxation Point
% NR = relaDura/sampleTime;

% index by time [second]
startTime = 0;                                                    % [second]
endTime = pusleStart+ (pulseDura+relaDura)*2*numel(PulseStrenth); % [second]
t = startTime:sampleTime:endTime;


%Voc = ones( TP, 1 )*2.2;           % oen-circuit voltage
load('Vinitial.mat');
Vini = VSTART(1,1);                        % open-circuit voltage

Vbat = zeros( NP+NR+NP+NR ,1 );     % terminal voltage
I    = zeros( NP+NR+NP+NR , 1 );

Voc = zeros( NP, 1 );
V0 = zeros( NP, 1 );                % Pulse voltage of R0
V1 = zeros( NP, 1 );                % Pulse voltage of 1st RC-pair 
V2 = zeros( NP, 1 );                % Pulse voltage of 2nd RC-pair

VocR = zeros( NR, 1 );
V0R = zeros( NR, 1 );               % Relaxation voltage of R0 
V1R = zeros( NR, 1 );               % Relaxation voltage of 1st RC-pair 
V2R = zeros( NR, 1 );               % Relaxation voltage of 2nd RC-pair 

 Vbat (startTime+1:StartID,1) = Vini.*ones(StartID,1)  ;    %Voc(startTime+1:StartID,1)
 I (startTime+1:StartID,1) = zeros(StartID,1)  ;            %Voc(startTime+1:StartID,1)
 Vhistory = [ ];
 Ihistory = [ ];
 
 %% EHPPC Start
 for numPulse = 1 : numel(PulseStrenth)
     
     pulseI = PulseStrenth(numPulse);
     
  
    % disch-Pulse Start
    for i = 1 : NP
            Voc (i,1)  =   Vini  - (PP(i)*pulseI)/C0;
            V0  (i)  =   R0*pulseI;
            V1  (i)  =   R1.*pulseI *(1- exp(-PP(i)/(R1*C1)));
            V2  (i)  =   R2.*pulseI *(1- exp(-PP(i)/(R2*C2)));    
    end
    
    Vbat ( StartID+1 : StartID+NP, 1 ) =  Voc - V0 - V1- V2;  %Voc(StartID+1 :StartID+NP,1)
    I    ( StartID+1 : StartID+NP ,1 ) = -pulseI ;
    
    % Relaxation      
    for i = 1 : NR    
            VocR (i,1)  =   Voc(end);
            V0R (i)  =    R0*0;
            V1R(i)   =    V1 (end).* exp(-RP(i)/(R1*C1));
            V2R(i)   =    V2 (end).* exp(-RP(i)/(R2*C2));
    end
    
    Vbat ( StartID+NP+1 : StartID+NP+NR, 1  ) = VocR - V0R - V1R- V2R; % Voc( StartID+NP+1 : StartID+NP+NR, 1)
     I   ( StartID+NP+1 : StartID+NP+NR, 1  ) = 0;


 
    % ch-Pulse Start 
    for i = 1 : NP
            Voc (i)  =   Voc(end) - (PP(i)*(-pulseI))/C0;
            V0 (i)   =   R0*(-pulseI);
            V1 (i)   =   R1.*(-pulseI) *(1- exp(-PP(i)/(R1*C1)));
            V2 (i)   =   R2.*(-pulseI) *(1- exp(-PP(i)/(R2*C2)));   
      
    end  
    
    Vbat  (StartID+NP+NR+1:StartID+NP+NR+NP,1) =  Voc - V0- V1- V2; %Voc(StartID+NP+NR+1 :StartID+NP+NR+NP,1)
     I    ( StartID+NP+NR+1:StartID+NP+NR+NP,1) = pulseI ; 

    
    % Relaxation      
    for i = 1 : NR 
        
            VocR (i)  =   Voc(end);
            V0R (i) =   R0*0;
            V1R(i) =    V1 (end).* exp(-RP(i)/(R1*C1));
            V2R(i) =    V2 (end).* exp(-RP(i)/(R2*C2));
    end
  
    Vbat (StartID+NP+NR+NP+1 :StartID+NP+NR+NP+NR,1)  = VocR - V0R - V1R- V2R; % Voc(StartID+NP+NR+NP+1 :StartID+NP+NR+NP+NR,1)
    I   ( StartID+NP+NR+NP+1 :StartID+NP+NR+NP+NR,1 ) = 0;
    
    Vhistory = [Vhistory Vbat'];
    Ihistory = [Ihistory I'];
    
   %Voc = Vbat(end);  % open-circuit voltage
    
   Voc = Vbat(end);
   Vbat = zeros( NP+NR+NP+NR ,1 );  
   I = zeros( NP+NR+NP+NR ,1 );    % terminal voltage
   StartID = 0;
   % Vbat (StartID+NP+NR+NP+NR : end ,1) =  Voc (StartID+NP+NR+NP+NR : end ,1);
 end



    Vbat =   transpose(Vhistory);
    %Vbat =  Vbat (1:63001);
    %Vbat = roundn( Vhistory, -5 );
    %Vbat0 = roundn( Vbat0, -5 );
    %Vbat0 =Vbat0(630000:end);
    j = 1;
  rmse = sqrt((sum((Vbat-Vbat0).^2))/size(Vbat0,1))
  error =  Vbat0-  Vbat;
  [maxError, maxIndex] = max(error)
    
figure;
plot(error);
xticks(startTime/sampleTime:1/sampleTime*60*10:endTime/sampleTime)
xticklabels({' '})
xlabel('Time(10 minute)');
ylabel('error(Vbat0-Vbat)');

figure;
sc = plot(Vbat0,'r');
% datatip(sc, 739500, 3.7408 );
% datatip(sc, 667500,  3.7101 );
hold on;
sg = plot(Vhistory,'b');
% datatip(sg,667500,  3.6955 );
hold on; 
xticks(startTime/sampleTime:1/sampleTime*60*10:endTime/sampleTime)
xticklabels({' '})
xlabel('Time(10 minute)');
legend('Vbat0', 'Vbat');
ylabel("V")

end
