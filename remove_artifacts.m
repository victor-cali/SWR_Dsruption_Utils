function [HPC,PFC,HPC_filtered,PFC_filtered]=remove_artifacts(HPC,PFC,fn,tr)

%Bandpass the signals before removing artifacts.
Wn1=[100/(fn/2) 300/(fn/2)]; % Cutoff=100-300 Hz
[b1,a1] = butter(3,Wn1,'bandpass'); 
HPC_filtered=filtfilt(b1,a1,HPC);
PFC_filtered=filtfilt(b1,a1,PFC);

TOT = abs(HPC) + abs(PFC); % Sum of amplitudes ==> To visually assess artifacts, as they will appear in every channel and add up

L = length(HPC(:,1));
% 
% figure
% plot(linspace(duration([0 0 0]),duration([0 0 L/fn]),L), TOT);

% Artifacts
% tr = 5670; %Visual threshold. Adjust it after viewing the signal.
outliers = false(L,1);
index = 1;
while index<L
    if TOT(index)>=tr
        if (index-300)>=0
        outliers((index-300):index+1999) = ones(2300,1); %When we detect an artifact, remove 300 datapoints prior (build up of artifact) and 2 sec of datapoints after.
        else
        outliers(1:index+1999) = ones(length(1:index+1999),1); %When we detect an artifact, remove 300 datapoints prior (build up of artifact) and 2 sec of datapoints after.
        end
        index = index+2000;
    else
        index = index +1;
    end
end


%Filter out artifacts & replace with the mean of the channels' medians

HPC(outliers,:) = mean(median(HPC)); 
PFC(outliers,:) = mean(median(PFC));

HPC_filtered(outliers,:) = mean(median(HPC_filtered));
PFC_filtered(outliers,:) = mean(median(PFC_filtered));
end
