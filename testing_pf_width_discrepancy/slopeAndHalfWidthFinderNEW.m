function [peak,positPeak,TopW,HW,Slopes1090,INFR,OUTFR]=slopeAndHalfWidthFinderNEW(RF,SmoothRF,ResampleRate)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    
    for cell=1:size(SmoothRF,1)
        if(0<sum(SmoothRF(cell,:)))
        [I,J]=max(SmoothRF(cell,:),[],2);        
        positPeak(cell)=J;
        peak(cell)=I;
        if((positPeak(cell)<size(SmoothRF,2))&&(1<positPeak(cell)))
            IndHW1(cell)=find(abs(SmoothRF(cell,[1:positPeak(cell)-1])-peak(cell)*0.5)==min(abs(SmoothRF(cell,[1:positPeak(cell)-1])-peak(cell)*0.5)));
            IndHW2(cell)=positPeak(cell)+find(abs(SmoothRF(cell,[positPeak(cell)+1:end])-peak(cell)*0.5)==min(abs(SmoothRF(cell,[positPeak(cell)+1:end])-peak(cell)*0.5)));
            IndTOP1(cell)=find(abs(SmoothRF(cell,[1:positPeak(cell)-1])-peak(cell)*0.95)==min(abs(SmoothRF(cell,[1:positPeak(cell)-1])-peak(cell)*0.95)));
            IndTOP2(cell)=positPeak(cell)+find(abs(SmoothRF(cell,[positPeak(cell)+1:end])-peak(cell)*0.95)==min(abs(SmoothRF(cell,[positPeak(cell)+1:end])-peak(cell)*0.95)));
            Ind101(cell)=find(abs(SmoothRF(cell,[1:positPeak(cell)-1])-peak(cell)*0.1)==min(abs(SmoothRF(cell,[1:positPeak(cell)-1])-peak(cell)*0.1)));
            Ind102(cell)=positPeak(cell)+find(abs(SmoothRF(cell,[positPeak(cell)+1:end])-peak(cell)*0.1)==min(abs(SmoothRF(cell,[positPeak(cell)+1:end])-peak(cell)*0.1)));
            Ind901(cell)=find(abs(SmoothRF(cell,[1:positPeak(cell)-1])-peak(cell)*0.9)==min(abs(SmoothRF(cell,[1:positPeak(cell)-1])-peak(cell)*0.9)));
            Ind902(cell)=positPeak(cell)+find(abs(SmoothRF(cell,[positPeak(cell)+1:end])-peak(cell)*0.9)==min(abs(SmoothRF(cell,[positPeak(cell)+1:end])-peak(cell)*0.9)));
        end
        if(positPeak(cell)==1)
            IndHW2(cell)=1+find(abs(SmoothRF(cell,[2:end])-peak(cell)*0.5)==min(abs(SmoothRF(cell,[2:end])-peak(cell)*0.5)));
            IndHW1(cell)=nan;
            IndTOP2(cell)=1+find(abs(SmoothRF(cell,[2:end])-peak(cell)*0.95)==min(abs(SmoothRF(cell,[2:end])-peak(cell)*0.95)));
            IndTOP1(cell)=nan;
            Ind102(cell)=1+find(abs(SmoothRF(cell,[2:end])-peak(cell)*0.1)==min(abs(SmoothRF(cell,[2:end])-peak(cell)*0.1)));
            Ind101(cell)=nan;
            Ind902(cell)=1+find(abs(SmoothRF(cell,[2:end])-peak(cell)*0.9)==min(abs(SmoothRF(cell,[2:end])-peak(cell)*0.9)));
            Ind901(cell)=nan;
        end
        if(positPeak(cell)==size(SmoothRF,2))
            IndHW1(cell)=find(abs(SmoothRF(cell,[1:end-1])-peak(cell)*0.5)==min(abs(SmoothRF(cell,[1:end-1])-peak(cell)*0.5)));
            IndHW2(cell)=nan;
            IndTOP1(cell)=find(abs(SmoothRF(cell,[1:end-1])-peak(cell)*0.95)==min(abs(SmoothRF(cell,[1:end-1])-peak(cell)*0.95)));
            IndTOP2(cell)=nan;
            Ind101(cell)=find(abs(SmoothRF(cell,[1:end-1])-peak(cell)*0.1)==min(abs(SmoothRF(cell,[1:end-1])-peak(cell)*0.1)));
            Ind102(cell)=nan;
            Ind901(cell)=find(abs(SmoothRF(cell,[1:end-1])-peak(cell)*0.9)==min(abs(SmoothRF(cell,[1:end-1])-peak(cell)*0.9)));
            Ind902(cell)=nan;
        end
        
        %Finding the the width of the top 95% of the receptive field
%         MinTop=[];MaxTop=[];
%         if(0.2*threshold(cell)<(peak(cell)-vall(cell)))
            if(isnan(IndTOP1(cell))&&~isnan(IndTOP2(cell)))
               TopW(cell)=2*(IndTOP2(cell)-positPeak(cell))/ResampleRate+rand()*0.1; 
            elseif(isnan(IndTOP2(cell))&&~isnan(IndTOP1(cell)))
               TopW(cell)=2*(positPeak(cell)-IndTOP1(cell))/ResampleRate+rand()*0.1; 
            elseif(~isnan(IndTOP2(cell))&&~isnan(IndTOP1(cell)))
               TopW(cell)=(IndTOP2(cell)-IndTOP1(cell))/ResampleRate+rand()*0.1;
            else   
               TopW(cell)=nan;
            end
            %Finding the the halfwidth of the receptive field
            if(isnan(IndHW1(cell))&&~isnan(IndHW2(cell)))
               HW(cell)=2*(IndHW2(cell)-positPeak(cell))/ResampleRate;
               INFR(cell)=nanmean(SmoothRF(cell,[positPeak(cell):IndHW2(cell)]));
               OUTFR(cell)=nanmean(SmoothRF(cell,[IndHW2(cell)+1:end]));
            elseif(isnan(IndHW2(cell))&&~isnan(IndHW1(cell)))
               HW(cell)=2*(positPeak(cell)-IndHW1(cell))/ResampleRate;
               INFR(cell)=nanmean(SmoothRF(cell,[IndHW1(cell):positPeak(cell)]));
               OUTFR(cell)=nanmean(SmoothRF(cell,[1:IndHW1(cell)-1]));
            elseif(~isnan(IndHW2(cell))&&~isnan(IndHW1(cell)))
               HW(cell)=(IndHW2(cell)-IndHW1(cell))/ResampleRate;
               INFR(cell)=nanmean(SmoothRF(cell,[IndHW1(cell):IndHW2(cell)]));
               OUTFR(cell)=nanmean([SmoothRF(cell,[1:IndHW1(cell)-1]),SmoothRF(cell,[IndHW2(cell)+1:end])]);
            else   
               HW(cell)=nan;
               INFR(cell)=nan;
               OUTFR(cell)=nan;
            end
            %Finding the 10% - 90% slope
            if(isnan(Ind901(cell))&&~isnan(Ind101(cell)))
               Slopes10901(cell)=(2*(Ind101(cell)-positPeak(cell))/(peak(cell)-SmoothRF(cell,Ind101(cell))))/ResampleRate; 
            elseif(~isnan(Ind901(cell))&&~isnan(Ind101(cell)))
               Slopes10901(cell)=((Ind901(cell)-Ind101(cell))/(SmoothRF(cell,Ind901(cell))-SmoothRF(cell,Ind101(cell))))/ResampleRate;
            else   
               Slopes10901(cell)=nan;
            end
            if(isnan(Ind902(cell))&&~isnan(Ind102(cell)))
               Slopes10902(cell)=(2*(Ind102(cell)-positPeak(cell))/(peak(cell)-SmoothRF(cell,Ind102(cell))))/ResampleRate; 
            elseif(~isnan(Ind902(cell))&&~isnan(Ind102(cell)))
               Slopes10902(cell)=((Ind902(cell)-Ind102(cell))/(-SmoothRF(cell,Ind902(cell))+SmoothRF(cell,Ind102(cell))))/ResampleRate;
            else   
               Slopes10902(cell)=nan;
            end
            Slopes1090(cell)=nanmean([Slopes10901(cell),Slopes10902(cell)]);
        else
            peak(cell)=nan;
            positPeak(cell)=nan;
            TopW(cell)=nan;
            HW(cell)=nan;
            Slopes1090(cell)=nan;
            INFR(cell)=nan;
            OUTFR(cell)=nan;
        end
            
    end
    positPeak=positPeak/ResampleRate;
    
    

end

