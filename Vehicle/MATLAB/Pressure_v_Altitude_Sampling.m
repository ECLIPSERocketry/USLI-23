function pv = Pressure_v_Altitude_Sampling(Zbo)

% Altitude and Pressure Tables from 

%John D Anderson, Fundamentals of Aerodynamics, 5th Edition
%ISBN 978-0-07-339810-5

alt= [0,500,1000,1500,2000,2500,3000,3499,3999,4499,4999,5499,5998,6498];
P_Table = 1000*[2.1162,2.0783,2.0409,2.0040,1.9677,1.9319,...
    1.8967,1.8619,1.8277,1.7941,1.7609,1.7282,1.6960,1.6643];


for i = [1:length(Zbo)]

    B = [find((alt < Zbo(i)),1,"last"),find((alt > Zbo(i)),1)];

    if (B(2)-B(1)) ~= 1
        pv(i) = P_Table(mean(B));
    else
        pv(i) = P_Table(B(1)) + (P_Table(B(2))-P_Table(B(1)))*(Zbo(i)-alt(B(1)))/(alt(B(2))-alt(B(1)));
    end
end