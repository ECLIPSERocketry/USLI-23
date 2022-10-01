function rho = Density_v_Altitude_Sampling(Zbo)

% Altitude and Density Tables from 

%John D Anderson, Fundamentals of Aerodynamics, 5th Edition
%ISBN 978-0-07-339810-5

alt= [0,500,1000,1500,2000,2500,3000,3499,3999,4499,4999,5499,5998,6498];
rho_Table = 0.001*[2.3769,2.3423,2.3081,2.2743,2.2409,2.2079,2.1752,...
    2.1429,2.1110,2.0794,2.0482,2.0174,1.9869,1.9567];


for i = [1:length(Zbo)]

    B = [find((alt < Zbo(i)),1,"last"),find((alt > Zbo(i)),1)];

    if (B(2)-B(1)) ~= 1
        rho(i) = rho_Table(mean(B));
    else
        rho(i) = rho_Table(B(1)) + (rho_Table(B(2))-rho_Table(B(1)))*(Zbo(i)-alt(B(1)))/(alt(B(2))-alt(B(1)));
    end
end