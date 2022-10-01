function DrF = Nakkas_DragReduction_Sampling_fz(N)

DiN = [0:20:980];
DrF_Table = 0.001*[1000,955,915,880,850,817,790,760,735,718,700,680,665,...
                    645,632,617,605,595,580,572,562,555,545,540,530,520,...
                    512,500,497,485,480,471,465,458,450,440,435,425,420,...
                    415,410,405,400,397,394,394,394,394,394,394];

% Copy paste lines 3-17 and uncomment for visual proof of sampling.
% Compare to graph on page 4 of https://www.nakka-rocketry.net/articles/altcalc.pdf

% plot(DiN,DrF_Table,'-*k')
% hold 'on'
% plot(DiN,DrF_Table,'b')
% xlim([0 1000])
% ylim([0 1.1])
% hold 'off'

for i = [1:length(N)]

    if N(i) > 979
        DrF(i) = 0.394;
    else
        B = [find((DiN < N(i)),1,"last"),find((DiN > N(i)),1)];

        if (B(2)-B(1)) ~= 1
            DrF(i) = DrF_Table(mean(B));
        else
            DrF(i) = DrF_Table(B(1)) + (DrF_Table(B(2))-DrF_Table(B(1)))*(N(i)-DiN(B(1)))/(DiN(B(2))-DiN(B(1)));
        end
    end
end