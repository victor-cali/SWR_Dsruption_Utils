
% %Load ADRITOOLS
% addpath(genpath('C:\Users\students\Documents\Swatantra\ADRITOOLS-master'))
% 
% %Load CorticoHippocampal
% addpath(genpath('C:\Users\students\Documents\Swatantra\CorticoHippocampal'))

% Path with the rat's data.
ratpath='/mnt/genzel/Rat/SWRDisruptionPlusMaze/plusmaze_toilet_data_correct/rat6';
%Load RGS14 github
addpath(genpath('/home/adrian/Downloads/LFP_RGS14'))
addpath('/home/adrian/Desktop/SWR_Dsruption_Utils')

ss=3;
fn=1000;
power_line_noise = 50 : 50 : 300;
cd(ratpath)

G=getfolder();

for j=1:length(G) % Iterate across study days.
cd(ratpath)
cd(G{j})

prepost=getfolder();
     for i=1: length(prepost)

    cd(ratpath)
    cd(G{j})

        cd(prepost{i})

        %Read brain areas and load states file.
    %% HPC     
    HPC=dir(strcat('*','HPC','*.mat'));
    HPC=HPC.name;
    HPC=load(HPC);
    HPC=getfield(HPC,'HPC_ripple');
    HPC=HPC.*(0.195);
    HPC = notch_filter2(HPC, fn, power_line_noise);

    Cortex=dir(strcat('*','PFC','*.mat'));
    Cortex=Cortex.name;
    Cortex=load(Cortex);
    Cortex=getfield(Cortex,'PFC');
    Cortex=Cortex.*(0.195);
    Cortex = notch_filter2(Cortex, fn, power_line_noise);

    %Load sleep scoring
    A = dir('*states*.mat');
    A={A.name};

    if  ~isempty(A)
           cellfun(@load,A);
    else
          error('No Scoring found')    
    end


    [sd_swr,signal2_hpc,signal2_pfc]=find_std(HPC,Cortex,states,ss,fn);
    Sd_Swr.sd5_hpc_co(j,i)=sd_swr.sd5_hpc_co;
    Sd_Swr.sd5_pfc_co(j,i)=sd_swr.sd5_pfc_co;

    %Concatenate NREM bouts
    co_pfc=cell2mat(signal2_pfc);
    co_hpc=cell2mat(signal2_hpc);
    TOT{i} = abs(co_pfc) + abs(co_hpc); % Sum of amplitudes ==> To visually assess artifacts, as they will appear in every channel and add up
    
    end
 if length(prepost)==1 %If there is only one session.
     Sd_Swr.sd5_hpc_co(j,i+1)=NaN;
     Sd_Swr.sd5_pfc_co(j,i+1)=NaN;
 end
% L = length([TOT{1}.' TOT{2}.']);
% close all
% allscreen()
% plot(linspace(duration([0 0 0]),duration([0 0 L/fn]),L), [TOT{1}.' TOT{2}.'])
% xlabel('Time')
% ylabel('Select threshold for artifacts')
% title('Concatenated NREM of pre and post')

artifact_th(j)=10*(std([TOT{1}.' TOT{2}.']))+mean([TOT{1}.' TOT{2}.'])

end
thresholds_perday_hpc=nanmean(Sd_Swr.sd5_hpc_co,2);
thresholds_perday_cortex=nanmean(Sd_Swr.sd5_pfc_co,2);

tr(1)=mean(thresholds_perday_hpc);
tr(2)=mean(thresholds_perday_cortex);

offset1={'0'};
offset2={'0'};

D1 = round(tr(1) + str2num(cell2mat(offset1)));
D2 = round(tr(2) + str2num(cell2mat(offset2)));

%% Detect and store all events.

All_events=[];
All_str=[];
for j=1:length(G) % Iterate across study days.
cd(ratpath)
cd(G{j})

prepost=getfolder();
     for i=1: length(prepost)

    cd(ratpath)
    cd(G{j})

        cd(prepost{i})

        %Read brain areas and load states file.
    %% HPC     
    HPC=dir(strcat('*','HPC','*.mat'));
    HPC=HPC.name;
    HPC=load(HPC);
    HPC=getfield(HPC,'HPC_ripple');
    HPC=HPC.*(0.195);
    HPC = notch_filter2(HPC, fn, power_line_noise);

    Cortex=dir(strcat('*','PFC','*.mat'));
    Cortex=Cortex.name;
    Cortex=load(Cortex);
    Cortex=getfield(Cortex,'PFC');
    Cortex=Cortex.*(0.195);
    Cortex = notch_filter2(Cortex, fn, power_line_noise);

    %Load sleep scoring
    A = dir('*states*.mat');
    A={A.name};

    if  ~isempty(A)
           cellfun(@load,A);
    else
          error('No Scoring found')    
    end
    
HPC_backup=HPC; %This is just for visualizing artifacts
Cortex_backup=Cortex;    
% Artifact removal
[HPC,Cortex,HPC_filtered,PFC_filtered]=remove_artifacts(HPC,Cortex,fn,artifact_th(j));

% [swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc,~,sig_cortex] = swr_check_thr(HPC,Cortex,states,ss,D1,D2,fn);
[swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc,~,sig_cortex] = detect_ripples(HPC_filtered,PFC_filtered,states,ss,D1,D2,fn);


['Found ' num2str(sum(s_hpc)) ' hippocampal ripples']
['Found ' num2str(sum(s_pfc)) ' cortical ripples']


%Get waveforms of cortical ripples and store them
waveforms=sig_cortex(~cellfun('isempty',sig_cortex));
waveforms=[waveforms{:}];
All_events=[All_events waveforms];

if i==1
    All_str.([G{j} '_presleep' ])=waveforms;   
else
    All_str.([G{j} '_post' ])=waveforms;       
end

     end
end

%% Generate table with counts and stores it.

fields = fieldnames(All_str)
Counts=[];
for field_id=1:numel(fields)
    si=All_str.(fields{field_id});
    [sa_mixed,si_mixed,th]=freq_specs(si,fn);
    close all
    Counts(field_id,:)=[length(si_mixed.g1) length(si_mixed.g2) length(sa_mixed.g1) length(sa_mixed.g2)];
end

T=table;
T.Variables=    [fields num2cell(Counts)];
T.Properties.VariableNames=[{'Trial'};{'Instant Slow'};{'Instant Fast'};{'Meanfreq Slow'};{'Meanfreq Fast'};];    
xo
cd(ratpath)
writetable(T,strcat('Rat4_slow_fast_counts.xls'),'Sheet',1,'Range','A2:Z50')  
save('Rat4_All_events.mat','All_events')

%% Plot distribution with all events and store thresholds.
si=[All_events]; %take all events.

[sa_mixed,si_mixed,th]=freq_specs(si,fn);

printing_image('Rat4_cortical_ripples')
close all
save('Rat4_TH.mat','th')

%% THIS IS THE END OF THE DETECTION
%The NEXT PART IS FOR PLOTTING PURPOSES.
    

%% Find epoch with most detections
max_length=cellfun(@length,swr_pfc(:,1));
N=max_length==max(max_length);


hpc=V_hpc{N};
pfc=V_pfc{N};
hpc2=signal2_hpc{N};
pfc2=signal2_pfc{N};
n=find(N);
if length(n)>1
    'Multiple epochs with same number of events'
end
n=n(1);
%%
prompt = {'Select window length (msec):','Brain area:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1000','PFC'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
 win_len=str2num(answer{1});
 BR=answer{2};
 %%

close all
plot((1:length(hpc))./1000,5.*zscore(hpc)+100,'Color','black')
hold on
plot((1:length(pfc))./1000,5.*zscore(pfc)+150,'Color','black')
xlabel('Time (Seconds)')


plot((1:length(hpc2))./1000,5.*zscore(hpc2)+220,'Color','black')
plot((1:length(pfc2))./1000,5.*zscore(pfc2)+290,'Color','black')


yticks([100 150 220 290])
yticklabels({'HPC','PFC','HPC (Bandpassed)',['PFC' '(Bandpassed)']})
b=gca;
b.FontSize=12;


if strcmp(BR,'PFC')
    sn=swr_pfc{n,3};
else
    sn=swr_hpc{n,3};
end

if isempty(sn)
    errordlg('No Events found','Error');
    xo
end


prompt = {['Select event ID number. Max value:' num2str(length(sn))]};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1'};
answer2 = inputdlg(prompt,dlgtitle,dims,definput);
answer2=str2num(answer2{1});

 
if ~ isempty(swr_hpc{n})
    stem([swr_hpc{n,3}],ones(length([swr_hpc{n}]),1).*250,'Color','blue') %(HPC)
end
hold on
if ~ isempty(swr_pfc{n})
    stem([swr_pfc{n,3}],ones(length([swr_pfc{n}]),1).*250,'Color','red')%Seconds (Cortex)
end

% xlim([sn(answer2)-win_len/1000 sn(answer2)+win_len/1000])

