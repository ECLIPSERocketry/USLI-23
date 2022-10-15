% Master Calculator for UC Davis Eclipse Rocketry Club.--------------------
% For use in the 2022 USLI competition.------------------------------------
% By Alec R. Schrader------------------------------------------------------


% Updated for Rocket 3.7---------------------------------------------------


clear, clc
format compact
format longG


% Universal Constants------------------------------------------------------

g = 32.22;                                                                  % Surface Gravitational Acceleration from Table D/2 'Engineering Mechanics Dynamics 9th Edition' by Meriam et al. (ft/s^2)
p = 0.07350/g;                                                              % Air density from Table D/1 'Engineering Mechanics Dynamics 9th Edition' by Meriam et al. (slug/ft^3)
ws = [0;5;10;15;20];                                                        % Wind Speed (mph)
la = 2.236;                                                                 % Launch Angle (Degrees)


% Conversions

lbts = (0.45359/14.594);                                                    % Pound mass to slugs from Table D/5 'Engineering Mechanics Dynamics 9th Edition' by Meriam et al.
mtf = 1.466667;                                                             % Miles to feet from https://www.unitconverters.net/speed/mph-to-feet-per-second.htm
ozts = 1/514.7848;                                                          % Ounces to slugs


% Vehicle Variables--------------------------------------------------------

F = 181;                                                                    % Average Thrust (lbf)
t = 2.17;                                                                   % Motor Burn time (seconds)

Cd = [0.60,0.65,0.66,0.70];                                                 % Dimensionless drag coefficient matrix

mml = 54.4;                                                                 % Motor Launch mass (oz)
mme = 23.7;                                                                 % Motor Empty mass (oz)
mp = (mml - mme);                                                           % Propellant mass (oz)

mn = (16.7);                                                                % Nose cone mass (oz)
ma1 = (19.4 + 10.7 + 1.2 + 5.8 + 1.04 + 0.6 + 14.4);                        % Airframe 1 mass (oz)
msb = (6.5 + 0.3 + 5.8 + 1.04 + 1 + 5.3 + 6.5 + 0.2);                       % Switch band mass (oz)
ma2 = (89.8 + 0.4 + 0.9 + 6 + mme);                                         % Airframe 2 mass (oz)
ma1sb = ma1 + msb;                                                          % Mass of airframe 1, with switch band machined onto it, (oz)
totalMass = mn + ma1sb + ma2 + mp;                                          % Total mass at Launch (oz)
md = totalMass-mp;                                                          % Rocket Dead mass (oz)
m = 0.5*mp + md;                                                            % Rocket Average mass (oz)

Dn = 4;                                                                     % Rocket Diameter (inches)
CdpD = 1.5;                                                                 % Coefficient of drag of drogue parachute
Dd = 12;                                                                    % Diameter of drogue parachute (inches)
CdpM = 2.2;                                                                 % Coefficient of drag of main parachute
Dm = 60;                                                                    % Diameter of main parachute (inches)

th = 0.079;                                                                 % Wall Thickness (inches)
di = 3.843;                                                                 % Inner Diameter (inches)

Ln = 16.5;                                                                  % inches
d = 4;                                                                      % inches
df = 4;                                                                     % inches
dr = 4;                                                                     % inches
Lt = 0;
Xp = 0;
Cr = 9;                                                                     % inches
Ct = 3;                                                                     % inches
S = 4;                                                                      % inches
Lf = sqrt(2.5^2 + S^2);                                                     % inches
R = d/2;                                                                    % inches
Xr = 5.5;                                                                   % inches
Xb = 78.5;                                                                  % inches
Nfins = 4;


% Simulation 1 data--------------------------------------------------------

sim1_apogee = [4956; 4944; 4916; 4878; 4843];                               % From Openrocket, simulated apogee matrix

sim1_descent_0 = (91.9 - 16.6);                                             % From Openrocket, simulated descent time
sim1_descent_5 = (92 - 16.6);                                               % From Openrocket, simulated descent time
sim1_descent_10 = (92.4 - 16.5);                                            % From Openrocket, simulated descent time
sim1_descent_15 = (91.7 - 16.5);                                            % From Openrocket, simulated descent time
sim1_descent_20 = (92.6 - 16.4);                                            % From Openrocket, simulated descent time

sim1_ghV = [16.2; 16.1; 16; 16.4; 16];                                      % Simulated ground hit velocity


% Simulation 2 data--------------------------------------------------------

sim2_apogee = [4956; 4943; 4906; 4865; 4848];                               % From Openrocket, simulated apogee matrix

sim2_descent_0 = (91.9 - 16.6);                                             % From Openrocket, simulated descent time
sim2_descent_5 = (92.3 - 16.6);                                             % From Openrocket, simulated descent time
sim2_descent_10 = (93.3 - 16.5);                                            % From Openrocket, simulated descent time
sim2_descent_15 = (93.6 - 16.4);                                            % From Openrocket, simulated descent time
sim2_descent_20 = (91.7 - 16.4);                                            % From Openrocket, simulated descent time

sim2_ghV = [16.2; 16.3; 16; 16.2; 16.1];                                    % Simulated ground hit velocity


% Simulation 3 data--------------------------------------------------------

sim3_apogee = [4956; 4944; 4913; 4871; 4829];                               % From Openrocket, simulated apogee matrix

sim3_descent_0 = (91.9 - 16.6);                                             % From Openrocket, simulated descent time
sim3_descent_5 = (91.9 - 16.6);                                             % From Openrocket, simulated descent time
sim3_descent_10 = (91.6 - 16.6);                                            % From Openrocket, simulated descent time
sim3_descent_15 = (92 - 16.5);                                              % From Openrocket, simulated descent time
sim3_descent_20 = (92.2 - 16.3);                                            % From Openrocket, simulated descent time

sim3_ghV = [16.2; 15.9; 16; 16.1; 16];                                      % Simulated ground hit velocity


% END OF VARIABLE INPUTS---------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Early excel output of vehicle metrics in original units (oz/in/etc.)-----

filename = 'Rocket_Hand_Calc.xls';
                                                                            % Rounded to two decimal places
DnR = round(Dn*100)/100;
totalMassR = round(totalMass*100)/100;
mpR = round(mp*100)/100;
mmlR = round(mml*100)/100;
tR = round(t*100)/100;
FR = round(F*100)/100;
gR = round(g*100)/100;
laR = round(la*100)/100;

vehM = {'Rocket Diameter';'Wet Mass';'Propellent Mass';...                  % Vehicle metrics labeling
        'Motor Launch Mass';'Burn Time';'Average Thrust';...
        'Gravitational Acceleration';'Launch Angle (From Vertical)'};
vehMnum = [DnR,totalMassR,mpR,mmlR,tR,FR,gR,laR];                           % Vehicle metrics numerical values
writecell({'Vehicle Metrics','Value'},filename,'Sheet',1,'Range','A1');
writecell(vehM,filename,'Sheet',1,'Range','A2');
writematrix(vehMnum.',filename,'Sheet',1,'Range','B2');


% Quick Calculations for later use-----------------------------------------

A = pi*(Dn/(2*12))^2;                                                       % Reference area of nose cone (ft^2)
Ad = pi*(Dd/(2*12))^2;                                                      % Reference area of drogue parachute (ft^2)
Am = pi*(Dm/(2*12))^2;                                                      % Reference area of main parachute (ft^2)

mdS= md*ozts;                                                               % Rocket Dead Mass (slug)
mS = m*ozts;                                                                % Rocket Average mass (slug)
mnS = mn*ozts;                                                              % Nose cone mass (slug)
ma1S = ma1*ozts;                                                            % Airframe 1 mass (slug)
msbS = msb*ozts;                                                            % Switch band (slug)
ma2S = ma2*ozts;                                                            % Airframe 2 mass (slug)
ma1sbS = ma1sb*ozts;                                                        % Mass of airframe 1, with switch band machined onto it (slug)
mSm = [mnS;ma1sbS;ma2S];                                                    % Mass Matrix (slugs)
msm = [mn;ma1sb;ma2];                                                       % Mass Matrix (oz)


% Barrowman Approximation for Center of Pressure---------------------------

Cnn = 2;                                                                    % Nose cone terms

Xn = 0.466*Ln;                                                              % Ogive -> Xn = 0.466*Ln;
                                                                            % Cone -> Xn = .666*Ln;

Cnt = 2*((dr/d)^2-(dr/d)^2);                                                % Conical Transition Term
Xt = Xp + Lt/3*(1+(1-df/dr)/(1-(df/dr)^2));

Cnf = (1+R/(S+R))*((4*Nfins*(S/d)^2)/(1+sqrt(1+(2*Lf/(Cr+Ct))^2)));         % Fin Terms
Xf = Xb + Xr/3*(Cr+2*Ct)/(Cr+Ct)+1/6*((Cr+Ct)-(Cr*Ct)/(Cr+Ct));

Cnr = Cnn + Cnt + Cnf;                                                      % Finding Center of Pressure

                                                                            % Xbar = (Cnn*Xn + Cnt*Xt+ Cnf*Xf)/Cnr
                                                                            % No transition

Xbar = (Cnn*Xn + Cnf*Xf)/Cnr;                                               % inches

Xbarsimulated = 69.177;                                                     % OpenRocket value in inches

percentDifference = (Xbar - Xbarsimulated)/((Xbar + Xbarsimulated))*200;


% Apogee Calculator Based on Nakka's Simplification------------------------

z1 = 0.5*((F/mS) - g)*t^2;                                                  % [ft]

v1 = sqrt(2*z1/mS*(F-mS*g));                                                % [ft/s]

Apogee_Vertical = F*z1/(mS*g);

N = (Cd*Dn^2*v1^2)/(24353*mdS);

fz = Nakkas_DragReduction_Sampling_fz(N);                                   % Calls the function that is an n=100 sampling of Nakka's chart

Apogee_Vertical = fz*Apogee_Vertical;                                       % Vertical Launch (ft)
Apogee_Actual = Apogee_Vertical*cosd(la);                                   % Accounting for launch angle (ft)


% Max V and Burnout Altitude from Nakka------------------------------------

fv = Nakkas_DragReduction_Sampling_fv(N);                                   % Calls the function that is an n=10 sampling of Nakka's chart
fzbo = Nakkas_DragReduction_Sampling_fzbo(N);                               % Calls the function that is an n=10 sampling of Nakka's chart

Zbo = fzbo * z1;                                                            % Motor Bunrout Altitude, from Nakka [ft]
vmax = fv * v1;                                                             % Max Velocity, from Nakka [ft/s]
pv = Pressure_v_Altitude_Sampling(Zbo);                                     % Atmospheric pressure at Burnout Altitude [lbf/ft^2]
rho = Density_v_Altitude_Sampling(Zbo);                                     % Atmospheric Density at Burnout Altitude [slugs/ft^3]


% Max Q and Stresses-------------------------------------------------------

Qmax = 0.5*rho.*vmax.^2;                                                    % Maximum Dynamic Pressure [lbf/ft^2]
ptotal = Qmax + pv;                                                         % Total Pressure at Max Q [lbf/ft^2]
Ohoop = ptotal*di/(2*th);                                                   % Maximum Hoop Stress [lbf/ft^2]
Olong = 2*Ohoop;                                                            % Maximum Longitudinal Stress [lbf/ft^2]


% Descent Time Calculator--------------------------------------------------

% Assuming the rocket is at rest at apogee, with the previously
% calculated Apogee_Actual

vD = sqrt((2*mdS*g)/(p*CdpD*Ad));                                           % Velocity with drogue parachute (ft/s)
vM = sqrt((2*mdS*g)/(p*CdpM*Am));                                           % Velocity with main parachute (ft/s)

dtD = (Apogee_Actual - 550)/vD;                                             % Descent time drogue (s)
dtM = 550/vM;                                                               % Descent time main (s)
dtActual = dtD + dtM;                                                       % Descent time actual (s)


% Drift Distance Calculator------------------------------------------------

drift = ws*mtf*dtActual;                                                    % Drift Matrix (ft)


% Kinetic Energy Calculator------------------------------------------------

ghvM = sqrt(vM^2 + (ws*mtf).^2);                                            % Ground hit velocity matrix
KEv = (0.5*mSm.*(vM)^2);                                                    % Kinetic energy vertical

KEh0 = 0.5*mSm*(ws(1)*mtf)^2;                                               % Kinetic energy horizontal, 0 mph
KEtotal0 = KEh0 + KEv;

KEh5 = 0.5*mSm*(ws(2)*mtf)^2;                                               % Kinetic energy horizontal, 5 mph
KEtotal5 = KEh5 + KEv;

KEh10 = 0.5*mSm*(ws(3)*mtf)^2;                                              % Kinetic energy horizontal, 10 mph
KEtotal10 = KEh10 + KEv;

KEh15 = 0.5*mSm*(ws(4)*mtf)^2;                                              % Kinetic energy horizontal, 15 mph
KEtotal15 = KEh15 + KEv;

KEh20 = 0.5*mSm*(ws(5)*mtf)^2;                                              % Kinetic energy horizontal, 20 mph
KEtotal20 = KEh20 + KEv;

KEtotal_matrix = [KEtotal0';KEtotal5';KEtotal10';KEtotal15';KEtotal20'];    % Kinetic energy matrix

% Simulation 1 Data calculations-------------------------------------------

sim1_descent_matrix = [sim1_descent_0;sim1_descent_5;sim1_descent_10;...    % Simulated descent time matrix
    sim1_descent_15;sim1_descent_20];
sim1_descent = [ws,sim1_descent_matrix];                                    % From Openrocket, simulated descent time matrix
sim1_drift = mtf*ws.*sim1_descent_matrix;
sim1_KE0 = 0.5*mSm'*sim1_ghV(1)^2;                                          % Simulated Kinetic Energy, 0 mph
sim1_KE5 = 0.5*mSm'*sim1_ghV(2)^2;                                          % Simulated Kinetic Energy, 5 mph
sim1_KE10 = 0.5*mSm'*sim1_ghV(3)^2;                                         % Simulated Kinetic Energy, 10 mph
sim1_KE15 = 0.5*mSm'*sim1_ghV(4)^2;                                         % Simulated Kinetic Energy, 15 mph
sim1_KE20 = 0.5*mSm'*sim1_ghV(5)^2;                                         % Simulated Kinetic Energy, 20 mph
sim1_KE_matrix = [sim1_KE0;sim1_KE5;sim1_KE10;sim1_KE15;sim1_KE20];         % Simulated Kinetic Energy Matrix


% Simulation 2 Data calculations-------------------------------------------

sim2_descent_matrix = [sim2_descent_0;sim2_descent_5;sim2_descent_10;...    % Simulated descent time matrix
    sim2_descent_15;sim2_descent_20];
sim2_descent = [ws,sim2_descent_matrix];                                    % From Openrocket, simulated descent time matrix
sim2_drift = mtf*ws.*sim2_descent_matrix;
sim2_KE0 = 0.5*mSm'*sim2_ghV(1)^2;                                          % Simulated Kinetic Energy, 0 mph
sim2_KE5 = 0.5*mSm'*sim2_ghV(2)^2;                                          % Simulated Kinetic Energy, 5 mph
sim2_KE10 = 0.5*mSm'*sim2_ghV(3)^2;                                         % Simulated Kinetic Energy, 10 mph
sim2_KE15 = 0.5*mSm'*sim2_ghV(4)^2;                                         % Simulated Kinetic Energy, 15 mph
sim2_KE20 = 0.5*mSm'*sim2_ghV(5)^2;                                         % Simulated Kinetic Energy, 20 mph
sim2_KE_matrix = [sim2_KE0;sim2_KE5;sim2_KE10;sim2_KE15;sim2_KE20];         % Simulated Kinetic Energy Matrix


% Simulation 3 Data calculations-------------------------------------------

sim3_descent_matrix = [sim3_descent_0;sim3_descent_5;sim3_descent_10;...    % Simulated descent time matrix
    sim3_descent_15;sim3_descent_20];
sim3_descent = [ws,sim3_descent_matrix];                                    % From Openrocket, simulated descent time matrix
sim3_drift = mtf*ws.*sim3_descent_matrix;
sim3_KE0 = 0.5*mSm'*sim3_ghV(1)^2;                                          % Simulated Kinetic Energy, 0 mph
sim3_KE5 = 0.5*mSm'*sim3_ghV(2)^2;                                          % Simulated Kinetic Energy, 5 mph
sim3_KE10 = 0.5*mSm'*sim3_ghV(3)^2;                                         % Simulated Kinetic Energy, 10 mph
sim3_KE15 = 0.5*mSm'*sim3_ghV(4)^2;                                         % Simulated Kinetic Energy, 15 mph
sim3_KE20 = 0.5*mSm'*sim3_ghV(5)^2;                                         % Simulated Kinetic Energy, 20 mph
sim3_KE_matrix = [sim3_KE0;sim3_KE5;sim3_KE10;sim3_KE15;sim3_KE20];         % Simulated Kinetic Energy Matrix


% Simulation Averages

simAp_matrix = [sim1_apogee,sim2_apogee,sim3_apogee]';                      % Simulated Apogee matrix
simAp_avg = [mean(simAp_matrix)]';                                          %#ok<NBRAK> % Simulated Apogee Average

simDes_matrix = [sim1_descent_matrix,sim2_descent_matrix,...                % Simulated Descent Matrix
    sim3_descent_matrix]';
simDes_avg = [mean(simDes_matrix)]';                                        %#ok<NBRAK> % Simulated Descent Average

simDrf_matrix = [sim1_drift,sim2_drift,sim3_drift]';                        % Simulated Drift Matrix
simDrf_avg = [mean(simDrf_matrix)]';                                        %#ok<NBRAK> % Simulated Drift Average


% Percent Difference-------------------------------------------------------

pd_ap = 200*abs(Apogee_Actual(3) - simAp_avg(1))./...                       % Percent Difference, Apogee
    (Apogee_Actual(3) + simAp_avg(1));
pd_dt = 200*abs(dtActual(3) - simDes_avg(1))./...                           % Percent Difference, Descent Time
    (dtActual(3) + simDes_avg(1));
pd_drf = 200*abs(drift(:,3) - simDrf_avg)./...                              % Percent Difference, Drift Distance
    (drift(:,3) + simDrf_avg);
pd_drf(1) = 0;

% Round to 2 Decimal Places------------------------------------------------

Apogee_Vertical=  round(Apogee_Vertical*100)/100;                           % Vertical Launch (ft)
Apogee_Actual = round(Apogee_Actual*100)/100;                               % Accounting for launch angle (ft)
dtActual = round(dtActual*100)/100;                                         % Descent time actual (s)
drift = round(drift*100)/100;                                               % Drift Matrix (ft)
vM = round(vM*100)/100;                                                     % Velocity with main parachute
KEtotal0 = round(KEtotal0*100)/100;                                         % Kinetic energy, 0 mph
KEtotal5 = round(KEtotal5*100)/100;                                         % Kinetic energy, 5 mph
KEtotal10 = round(KEtotal10*100)/100;                                       % Kinetic energy, 10 mph
KEtotal15 = round(KEtotal15*100)/100;                                       % Kinetic energy, 15 mph
KEtotal20 = round(KEtotal20*100)/100;                                       % Kinetic energy, 20 mph
ghvM = round(ghvM*100)/100;                                                 % Ground hit velocity main
sim1_descent = round(sim1_descent*100)/100;                                 % Simulation 1 descent time matrix
sim2_descent = round(sim2_descent*100)/100;                                 % Simulation 2 descent time matrix
sim3_descent = round(sim3_descent*100)/100;                                 % Simulation 3 descent time matrix
simAp_avg = round(simAp_avg*100)/100;                                       % Simulated Apogee Average
pd_ap = round(pd_ap*100)/100;                                               % Percent Difference, Apogee
pd_dt = round(pd_dt*100)/100;                                               % Percent Difference, Descent Time
pd_drf = round(pd_drf*100)/100;                                             % Percent Difference, Drift Distance


% Table Buillding----------------------------------------------------------

Tapogee = {'Cd';'Altitude [ft]'};                                           % Apogee title heading
Tdescent = {'Cd', 'Descent Time [s]'};                                      % Time of descent title heading
simAp = {'Wind Speed [mph]', 'Simulation 1 Altitude [ft]',...               % Simulated Apogee title heading
    'Simulation 2 Altitude [ft]','Simulation 3 Altitude [ft]'};
simDes = {'Wind Speed [mph]', 'Simulation 1 Descent Time [s]',...           % Simulated Descent title heading
    'Simulation 2 Descent Time [s]','Simulation 3 Descent Time [s]'};
simDrf = {'Wind Speed [mph]', 'Simulation 1 Drift Distance [ft]',...        % Simulated Drift Distance title heading
    'Simulation 2 Drift Distance [ft]','Simulation 3 Drift Distance [ft]'};
keTitle = {'Wind Speed [mph]','Ground Hit Velocity [ft/s]',...              % Kinetic Energy title heading
    'Nosecone Kinetic Energy [ft-lbf]',...
    'Airframe 1 Kinetic Energy [ft-lbf]',...
    'Airframe 2 Kinetic Energy [ft-lbf]'};
keSub = {'Nose cone';'Airframe 1';'Airframe 2'};                            % Kinetic energy subtitle
boTitle = {'Coefficient of Drag','Burnout Altitude','Burnout Velocity'};    % Burnout Title
qTitle = {'Coefficient of Drag','Maximum Dynamic Pressure [lbf/ft^2]',...   % Max Q and Stresses Title
    'Static Pressure [lbf/ft^2]','Hoop Stress [lbf/ft^2]',...
    'Longitudinal Stress [lbf/ft^2]'};


% Hand Calculated Excel Exporting------------------------------------------

writecell({'Hand Calculated';'Results'},filename,'Sheet',1,'Range','D1')
                                                                            % Vertical Apogee
writematrix('Vertical Altitude Calculations',filename,'Sheet',1,'Range' ...
    ,'E2');
writecell(Tapogee,filename,'Sheet',1,'Range','F3')
writematrix(Cd,filename,'Sheet',1,'Range','G3')
writematrix(Apogee_Vertical,filename,'Sheet',1,'Range','G4')
                                                                            % Actual Apogee
writematrix('Preliminary Apogee Results',filename,'Sheet',1,'Range','E6')  
writecell(Tapogee,filename,'Sheet',1,'Range','F7')
writematrix(Cd,filename,'Sheet',1,'Range','G7')
writematrix(Apogee_Actual,filename,'Sheet',1,'Range','G8')
                                                                            % Descent Time
writematrix('Hand Calculated Descent Time',filename,'Sheet',1,'Range', ...
    'E10')
writecell(Tdescent,filename,'Sheet',1,'Range','F11')
writematrix(Cd.',filename,'Sheet',1,'Range','F12')
writematrix(dtActual.',filename,'Sheet',1,'Range','G12')
                                                                            % Drift Distance
writematrix('Hand Calculated Drift Distance',filename,'Sheet',1,'Range' ...
    ,'E17')
writematrix('Wind Speed',filename,'Sheet',1,'Range','F19')
writematrix('Coefficient of Drag',filename,'Sheet',1,'Range','G18')
writematrix(ws,filename,'Sheet',1,'Range','F20')
writematrix(Cd,filename,'Sheet',1,'Range','G19')
writematrix(drift,filename,'Sheet',1,'Range','G20')
                                                                            % Kinetic Energy
writematrix('Kinetic Energy Calculation',filename,'Sheet',1,'Range','E26')

writecell({'Section','Mass [oz]'},filename,'Sheet',1,'Range','F27')
writecell(keSub,filename,'Sheet',1,'Range','F28')
writematrix(msm,filename,'Sheet',1,'Range','G28')

writecell(keTitle,filename,'Sheet',1,'Range','F32')
writematrix(ws,filename,'Sheet',1,'Range','F33')
writematrix(ghvM,filename,'Sheet',1,'Range','G33')
writematrix(KEtotal_matrix,filename,'Sheet',1,'Range','H33')
                                                                            % Burnout Alt and V
writecell(boTitle,filename,'Sheet',1,'Range','L7')
writematrix(Cd',filename,'Sheet',1,'Range','L8')
writematrix(Zbo',filename,'Sheet',1,'Range','M8')
writematrix(vmax',filename,'Sheet',1,'Range','N8')
                                                                            % Max Q and Stresses
writecell(qTitle,filename,'Sheet',1,'Range','L13')
writematrix(Cd',filename,'Sheet',1,'Range','L14')
writematrix(Qmax',filename,'Sheet',1,'Range','M14')
writematrix(pv',filename,'Sheet',1,'Range','N14')
writematrix(Ohoop',filename,'Sheet',1,'Range','O14')
writematrix(Olong',filename,'Sheet',1,'Range','P14')

writecell({'Center of Pressure','Length From Front End [in]'},...           % Center of Pressure
    filename,'Sheet',1,'Range','I26')
writematrix(Xbar,filename,'Sheet',1,'Range','J27')

% Simulation data excel exporting----------------------------------------

writematrix('Simulation Results',filename,'Sheet',2,'Range','A1')
                                                                            % Altitude
writematrix('Simulation Apogee',filename,'Sheet',2,'Range','B2')
writecell(simAp,filename,'Sheet',2,'Range','C3')
writematrix(ws,filename,'Sheet',2,'Range','C4')
writematrix(sim1_apogee,filename,'Sheet',2,'Range','D4')
writematrix(sim2_apogee,filename,'Sheet',2,'Range','E4')
writematrix(sim3_apogee,filename,'Sheet',2,'Range','F4')
writecell({'Average Simulated Altitude'},filename,'Sheet',2,'Range','G3')
writematrix(simAp_avg,filename,'Sheet',2,'Range','G4')
                                                                            % Descent Time
writematrix('Simulation Descent Time',filename,'Sheet',2,'Range','B10')
writecell(simDes,filename,'Sheet',2,'Range','C11')
writematrix(ws,filename,'Sheet',2,'Range','C12')
writematrix(sim1_descent_matrix,filename,'Sheet',2,'Range','D12')
writematrix(sim2_descent_matrix,filename,'Sheet',2,'Range','E12')
writematrix(sim3_descent_matrix,filename,'Sheet',2,'Range','F12')
writecell({'Average Descent Time [s]'},filename,'Sheet',2,'Range','G11')
writematrix(simDes_avg,filename,'Sheet',2,'Range','G12')
                                                                            % Drift Distance
writematrix('Simulation Drift Distance',filename,'Sheet',2,'Range','B18')
writecell(simDrf,filename,'Sheet',2,'Range','C19')
writematrix(ws,filename,'Sheet',2,'Range','C20')
writematrix(sim1_drift,filename,'Sheet',2,'Range','D20')
writematrix(sim2_drift,filename,'Sheet',2,'Range','E20')
writematrix(sim3_drift,filename,'Sheet',2,'Range','F20')
writecell({'Average Drift Distance [ft]'},filename,'Sheet',2,'Range','G19')
writematrix(simDrf_avg,filename,'Sheet',2,'Range','G20')
                                                                            % Kinetic Energy
writematrix('Simulation Kinetic Energy',filename,'Sheet',2,'Range','B26')

writecell({'Simulation Center of Pressure',...                              % Center of Pressure
    'Length From Front End [in]'},filename,'Sheet',2,'Range','F26')
writematrix(Xbarsimulated,filename,'Sheet',2,'Range','G27')

writecell({'Section','Mass [oz]'},filename,'Sheet',2,'Range','C27')
writecell(keSub,filename,'Sheet',2,'Range','C28')
writematrix(msm,filename,'Sheet',2,'Range','D28')

                                                                            % Simulation 1
writecell(keTitle,filename,'Sheet',2,'Range','C32')
writematrix(ws,filename,'Sheet',2,'Range','C33')
writematrix(sim1_ghV,filename,'Sheet',2,'Range','D33')
writematrix(sim1_KE_matrix,filename,'Sheet',2,'Range','E33')
                                                                            % Simulation 2
writecell(keTitle,filename,'Sheet',2,'Range','C39')
writematrix(ws,filename,'Sheet',2,'Range','C40')
writematrix(sim2_ghV,filename,'Sheet',2,'Range','D40')
writematrix(sim3_KE_matrix,filename,'Sheet',2,'Range','E40')
                                                                            % Simulation 3
writecell(keTitle,filename,'Sheet',2,'Range','C46')
writematrix(ws,filename,'Sheet',2,'Range','C47')
writematrix(sim1_ghV,filename,'Sheet',2,'Range','D47')
writematrix(sim1_KE_matrix,filename,'Sheet',2,'Range','E47')


% Percent Difference-------------------------------------------------------

writecell({'Percent Difference'},filename,'Sheet',3,'Range','A1')

writecell({'Apogee'},filename,'Sheet',3,'Range','B2')
writecell({'Percent Difference of Apogee'}, ...
    filename,'Sheet',3,'Range','C3')
writematrix(pd_ap,filename,'Sheet',3,'Range','C4')

writecell({'Descent Time'},filename,'Sheet',3,'Range','B6')
writecell({'Percent Difference of Descent Time'}, ...
    filename,'Sheet',3,'Range','C7')
writematrix(pd_dt,filename,'Sheet',3,'Range','C8')

writecell({'Drift Distance'},filename,'Sheet',3,'Range','B10')
writecell({'Wind Speed [mph]','Percent Difference of Drift Distance'}, ...
    filename,'Sheet',3,'Range','C11')
writematrix(ws,filename,'Sheet',3,'Range','C12')
writematrix(pd_drf,filename,'Sheet',3,'Range','D12')

writecell({'Center of Pressure'},filename,'Sheet',3,'Range','B18')
writecell({'Percent Difference of Center of Pressure'},filename,...
    'Sheet',3,'Range','C19')
writematrix(percentDifference,filename,'Sheet',3,'Range','C20')


% Nakka Output-------------------------------------------------------------

writecell({'Drag Influence Number'},filename,'Sheet',1,'Range','L1')
writecell({'fz'},filename,'Sheet',1,'Range','M1')
writecell({'fv'},filename,'Sheet',1,'Range','N1')
writecell({'fzbo'},filename,'Sheet',1,'Range','O1')
writematrix(N',filename,'Sheet',1,'Range','L2')
writematrix(fz',filename,'Sheet',1,'Range','M2')
writematrix(fv',filename,'Sheet',1,'Range','N2')
writematrix(fzbo',filename,'Sheet',1,'Range','O2')


%--------------------------------------------------------------------------
xlsAutoFitCol(filename,'Sheet3','A:Z')
xlsAutoFitCol(filename,'Sheet2','A:Z')
xlsAutoFitCol(filename,'Sheet1','A:Z')
winopen('Rocket_Hand_Calc.xls')
%--------------------------------------------------------------------------