function [ Vbat,Vini ] = batVoltResponse( P,VSTART )
%%
    global PULSE
    R0 = P(1);
    R1 = P(2);
    C1 = P(3);
    R2 = P(4);
    C2 = P(5);
    C0 = P(6);
    
%% EHPPC 參數設定5
    pusleStart = 0;     % [second]
    pulseDura  = 30.01;     % [second]
    relaDura   = 600;    % [second]
    sampleTime = 0.01;   % [second]

    PulseStrenth = PULSE;       % Strenth of the Pulses[A]

%% 初始化
%

    StartID = pusleStart/sampleTime;  % data point index
    PP = 0:sampleTime:pulseDura;    % Duration of pulse
    NP = pulseDura/sampleTime;     % Number of Pulse Point
    RP = 0:sampleTime:relaDura;     % Duration of relaxation
    NR = relaDura/sampleTime;      % Number of Relaxation Point
% NR = relaDura/sampleTime;
    
    startTime = 0;                      % [second]
    endTime = pusleStart+ (pulseDura+relaDura)*2*numel(PulseStrenth);      % [second]
    t = startTime:sampleTime:endTime;


    load('Vinitial.mat');
    Vini = VSTART(1,1);                       % open-circuit voltage
    
    
    V    = zeros( NP+NR+NP+NR , 1 );    % terminal voltage
    I    = zeros( NP+NR+NP+NR , 1 );
   
    Voc = zeros( NP, 1 );
    V0 = zeros( NP, 1 );    % Pulse voltage of R0
    V1 = zeros( NP, 1 );    % Pulse voltage of 1st RC-pair 
    V2 = zeros( NP, 1 );    % Pulse voltage of 2nd RC-pair
    
    VocR = zeros( NR, 1 );
    V0R = zeros( NR, 1 );   % Relaxation voltage of R0 
    V1R = zeros( NR, 1 );   % Relaxation voltage of 1st RC-pair 
    V2R = zeros( NR, 1 );   % Relaxation voltage of 2nd RC-pair 

     V (startTime+1:StartID,1) = Vini.*ones(StartID,1)  ; %Voc(startTime+1:StartID,1)
     I (startTime+1:StartID,1) = zeros(StartID,1)  ; %Voc(startTime+1:StartID,1)
     Vhistory = [ ];
     Ihistory = [ ];
 
      %% EHPPC Start
     for numPulse = 1 : numel(PulseStrenth)

         pulseI = PulseStrenth(numPulse);


        % disch-Pulse Start
        for i = 1 : NP
                Voc (i,1)   =   Vini  - (PP(i)*pulseI)/C0;
                V0 (i)      =   R0*pulseI;
                V1 (i)      =   R1.*pulseI *(1- exp(-PP(i)/(R1*C1)));
                V2 (i)      =   R2.*pulseI *(1- exp(-PP(i)/(R2*C2)));    
        end

        V ( StartID+1 : StartID+NP, 1 ) =  Voc - V0 - V1- V2; %Voc(StartID+1 :StartID+NP,1)
        I ( StartID+1 : StartID+NP ,1 ) = -pulseI ;

        % Relaxation      
        for i = 1 : NR      
                VocR (i,1)  =   Voc(end);
                V0R (i)     =   R0*0;
                V1R(i)      =   V1 (end).* exp(-RP(i)/(R1*C1));
                V2R(i)      =   V2 (end).* exp(-RP(i)/(R2*C2));
        end

        V ( StartID+NP+1 : StartID+NP+NR, 1  ) = VocR - V0R - V1R - V2R; %Voc( StartID+NP+1 : StartID+NP+NR, 1)
        I ( StartID+NP+1 : StartID+NP+NR, 1  ) = 0;



        % ch-Pulse Start 
        for i = 1 : NP
                Voc (i,1)   =   Voc(end) - (PP(i)*(-pulseI))/C0;
                V0 (i)      =   R0*(-pulseI);
                V1 (i)      =   R1.*(-pulseI) *(1- exp(-PP(i)/(R1*C1)));
                V2 (i)      =   R2.*(-pulseI) *(1- exp(-PP(i)/(R2*C2)));   

        end  

        V ( StartID+NP+NR+1:StartID+NP+NR+NP,1 ) =  Voc - V0- V1- V2; %Voc(StartID+NP+NR+1 :StartID+NP+NR+NP,1)
        I ( StartID+NP+NR+1:StartID+NP+NR+NP,1 ) = pulseI ; 


        % Relaxation      
        for i = 1 : NR      
                VocR (i,1)  =   Voc(end);
                V0R (i)     =   R0*0;
                V1R(i)      =    V1 (end).* exp(-RP(i)/(R1*C1));
                V2R(i)      =    V2 (end).* exp(-RP(i)/(R2*C2));
        end

        V ( StartID+NP+NR+NP+1 :StartID+NP+NR+NP+NR,1 )  = VocR - V0R - V1R - V2R; %Voc(StartID+NP+NR+NP+1 :StartID+NP+NR+NP+NR,1)
        I ( StartID+NP+NR+NP+1 :StartID+NP+NR+NP+NR,1 )  = 0;

        Vhistory = [Vhistory V'];
        Ihistory = [Ihistory I'];

       Voc = V(end);
       V = zeros( NP+NR+NP+NR ,1 );  
       I = zeros( NP+NR+NP+NR ,1 );    % terminal voltage
       StartID = 0;

     end
     
     
     Vbat = Vhistory;
      
end

