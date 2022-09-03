clear
clc
for k = 1:100   
tic
QQ = 50;     %SOC
sampleTime = 0.01;
FileName = ['C:\Users\lcy20\Desktop\DATA 23\SOC',num2str(QQ),'.xlsx'];
%%',num2str(QQ),'%%       
DataImport = readtable(FileName);
Data = table2array(DataImport);
[m, n] = size(Data);
%load Data
PlotFrom            = 2;                      % number of data point 
PlotTo              = m;                % number of data point
V = Data(PlotFrom:PlotTo, 1);
VSTART = Data(1, 1);        %SOC point's   NCV
I = Data(PlotFrom:PlotTo, 2);
stepTime =  Data(PlotFrom:PlotTo, 1);
save('Vmeasure.mat','V')       
save('Vinitial.mat','VSTART')
[ R0,R1,C1,R2,C2,C0,gBest,Fitness_gBest ] = pso( V );   %call pso.m
bestresult = gBest(100,1:6);
[ maxError,Vbat0, Vbat  ] = batVoltResNEW( R0,R1,C1,R2,C2,C0 ); %call batVoltResNEW
toc;
K(k,1) = R0;
K(k,2) = R1;
K(k,3) = C1;
K(k,4) = R2;
K(k,5) = C2;
K(k,6) = C0;
K(k,7) = Fitness_gBest(100,1);
K(k,8) = maxError;
K(k,9) = toc;
end