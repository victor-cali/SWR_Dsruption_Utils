function [HPC,PFC]=remove_artifacts(HPC,PFC,fn,tr)
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
        if 300>=index; low_lim = 1; else; low_lim = index-300; end
        if 1999>length(outliers); high_lim = length(outliers); else; high_lim = index+1999; end
        outliers(low_lim:high_lim) = ones(length(low_lim:high_lim),1); %When we detect an artifact, remove 300 datapoints prior (build up of artifact) and 2 sec of datapoints after. 
        index = index+2000;
    else
        index = index +1;
    end
end


%Filter out artifacts & replace with the mean of the channels' medians

HPC(outliers,:) = mean(median(HPC)); 
PFC(outliers,:) = mean(median(PFC));
end