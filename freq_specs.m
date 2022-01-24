function [si_mixed,sa_mixed,TH]=freq_specs(si,fn)
%Computes main features of events detected.
    if ~isempty(si)

        %Instantaneous frequency.
        x=cellfun(@(equis) mean(instfreq(equis,fn)) ,si,'UniformOutput',false);
        x=cell2mat(x);
        
%         th=gaussmix(x,Rat,tr);
        [th]=gaussmix_slow_fast(x);
        TH(1)=th;
        th=167;
        sa_mixed.g1=si(x<=th);
        sa_mixed.i1=find(x<=th);
        sa_mixed.g2=si(x>th);
        sa_mixed.i2=find(x>th);
            allscreen();
            subplot(2,2,1)
            histogram(x,[100:5:250]); title('Instantaneous Frequencies');xlabel('Frequency (Hz)');ylabel('Count')
%            histogram(x); title('Instantaneous Frequencies');xlabel('Frequency (Hz)');ylabel('Count')

            xline(th,'-',num2str(th),'LineWidth',1)            
            xlim([100 250])
%             ylim([0 500])

            %             ylim([0 35])
%             yticks([0:5:35])

%             histogram(x,[100:10:250],'Normalization','probability'); title('Instantaneous Frequencies');xlabel('Frequency (Hz)');ylabel('Probability')
%             ylim([0 0.5])
%             yticks([0:0.1:0.5])
        %Average frequency
        y=cellfun(@(equis) (meanfreq(equis,fn)) ,si,'UniformOutput',false);
        y=cell2mat(y);
%        th=gaussmix(y,Rat,tr);
         [th]=gaussmix_slow_fast(y);
        TH(2)=th;     
        th=155;
        si_mixed.g1=si(y<=th);
        si_mixed.i1=find(y<=th);
        si_mixed.g2=si(y>th);
        si_mixed.i2=find(y>th);
 
            subplot(2,2,2)
            histogram(y,[100:5:250]); title('Average Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
%            histogram(y); title('Average Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 

            xline(th,'-',num2str(th),'LineWidth',1)
            xlim([100 250])

%             ylim([0 50])
%             
%             yticks([0:10:50])
%             histogram(y,[100:10:250],'Normalization','probability'); title('Average Frequencies');xlabel('Frequency (Hz)');ylabel('Probability')
%             ylim([0 0.5]) 
%             yticks([0:0.1:0.5])            

%         yfreqmax=cellfun(@(equis) (freqmaxpeak(equis,fn)) ,si,'UniformOutput',false);
%         yfreqmax=cell2mat(yfreqmax);
%         
% %        th=gaussmix(yfreqmax,Rat,tr);
%          [th]=gaussmix_slow_fast(yfreqmax);
%                  TH(3)=th;
%                  th=155;
%         sp_mixed.g1=si(yfreqmax<=th);
%         sp_mixed.i1=find(yfreqmax<=th);
%         sp_mixed.g2=si(yfreqmax>th);
%         sp_mixed.i2=find(yfreqmax>th);

        
%              subplot(3,2,3)
% 
%             histogram(yfreqmax,[100:10:250]); title('Peak Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
% %            histogram(yfreqmax); title('Peak Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
% 
%             xline(th,'LineWidth',1)
%             xlim([100 250])
% %%
%         yfreqmax=cellfun(@(equis) (meanfreq2(equis,fn)) ,si,'UniformOutput',false);
%         yfreqmax=cell2mat(yfreqmax);
%         
% %        th=gaussmix(yfreqmax,Rat,tr);
%          [th]=gaussmix_slow_fast(yfreqmax);
%                  TH(4)=th;
%                  th=155;
%         sp2_mixed.g1=si(yfreqmax<=th);
%         sp2_mixed.i1=find(yfreqmax<=th);
%         sp2_mixed.g2=si(yfreqmax>th);
%         sp2_mixed.i2=find(yfreqmax>th);
% 
%              subplot(3,2,4)
% 
%             histogram(yfreqmax,[100:10:250]); title('Peak Frequencies. Method 2');xlabel('Frequency (Hz)');ylabel('Count') 
% %            histogram(yfreqmax); title('Peak Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
% 
%             xline(th,'LineWidth',1)
%             xlim([100 250])
%% sqrt
%         yfreqmax=cellfun(@(equis) (meanfreq2sqrt(equis,fn)) ,si,'UniformOutput',false);
%         yfreqmax=cell2mat(yfreqmax);
%         
% %        th=gaussmix(yfreqmax,Rat,tr);
%          [th]=gaussmix_slow_fast(yfreqmax);
%                  TH(5)=th;
%                  th=188;
%         spsqrt_mixed.g1=si(yfreqmax<=th);
%         spsqrt_mixed.i1=find(yfreqmax<=th);
%         spsqrt_mixed.g2=si(yfreqmax>th);
%         spsqrt_mixed.i2=find(yfreqmax>th);
% 
%              subplot(3,2,5)
% 
%             histogram(yfreqmax,[100:5:250]); title('Squared Power * Frequencies.');xlabel('Frequency (Hz)');ylabel('Count') 
% %            histogram(yfreqmax); title('Peak Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
% 
%             xline(th,'LineWidth',1)
%             xlim([100 250])


%%
            subplot(3,2,6)
            text(0.2,0.5,{['Inst:' num2str(length(sa_mixed.g1)) ' Slow HFOs'],...
                ['Inst:' num2str(length(sa_mixed.g2)) ' Fast HFOs'],...
                ['------'],...
                ['Average:' num2str(length(si_mixed.g1)) ' Slow HFOs'],...
                ['Average:' num2str(length(si_mixed.g2)) ' Fast HFOs'],...
%                 ['------'],...
%                 ['Peak:' num2str(length(sp_mixed.g1)) ' Slow HFOs'],...
%                 ['Peak:' num2str(length(sp_mixed.g2)) ' Fast HFOs'],...
%                 ['------'],...
%                 ['Peak2:' num2str(length(sp2_mixed.g1)) ' Slow HFOs'],...
%                 ['Peak2:' num2str(length(sp2_mixed.g2)) ' Fast HFOs'],...
%                 ['------'],...
%                 ['Sqrt:' num2str(length(spsqrt_mixed.g1)) ' Slow HFOs'],...
%                 ['Sqrt:' num2str(length(spsqrt_mixed.g2)) ' Fast HFOs'],...
%                 
                },'FontSize',14); axis off      
    else

        x=NaN;
        y=NaN;
        z=NaN;
        w=NaN;
        h=NaN;
        q=NaN;
        l=NaN;
        p=NaN;
        th=NaN;
        si_mixed.g1=NaN;
        si_mixed.i1=NaN;
        si_mixed.g2=NaN;
        si_mixed.i2=NaN;
        sa_mixed.g1=NaN;
        sa_mixed.i1=NaN;
        sa_mixed.g2=NaN;
        sa_mixed.i2=NaN;
        
        TH=NaN;
    end
    
%freqmaxpeak(x,fn)
end
