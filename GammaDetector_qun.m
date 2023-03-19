% function gammaevents
% based on cell data
% Output: GammaEvents
%        Column: start time, max, end time, frequency, cycles
function [GammaEvents] = GammaDetector_qun(C,savename,Fs,detchan,minuschan)
% nChannels is channel number in recordings
% detchan and minuschan is channel number in matlab,if in neuroscope-1
% if detchan is a vector, minuschan should be 
% minuschan(optinal) (physiological position) is far from detect channel for remove noise;
%parameter
% Fs = 500; %S/s
cutoffhighpass = 30;
cutofflowpass = 80; %LG 30-80;

window = 50; % ms
detectionlevel = 3; %SD, try different levels... erans paper 5,2
boundarylevel = 2; %SD, it does not make sense to put it too high if I'm going to use 5 / 7SD 5/7SD

%nChannels = size(C{1,1},2);
nepoch = length(C);
C1 = cell(nepoch,1);
%% pick up channel
if nargin <5
    for i = 1:nepoch
        C1{i} = C{i}(:,detchan);
    end
else
    if length(detchan)~=length(minuschan)
        error('channel length are not equal..');        
    else
        for i = 1:nepoch
            C1{i} = C{i}(:,detchan)-C{i}(:,minuschan);
        end
    end
end
clear C
%% detect gamma 30-80 Hz
gammaname = 'Gamma';
GammaEvents = cell(nepoch,1);
len = length(detchan);
for p = 1:nepoch
    C2 = C1{p};
    disp(['processing epoch ' num2str(p) '......'])
    %%           
    for ch = 1:len 
        %%
        % Filter channel
        fNorm = cutoffhighpass / (Fs /2); 
        [d,e] = butter(4, fNorm, 'high');
        temp= filtfilt(d, e, double(C2(:,ch)));  

        fNorm = cutofflowpass / (Fs /2);
        [d,e] = butter(4, fNorm, 'low');
        Filtchan= filtfilt(d, e, temp);   
        clear temp d e fNorm
        % calculate power in x ms windows
        power = zeros(size(Filtchan));
        winlen = size(Filtchan,1)- ceil(window*Fs/1000);
        for i = 1: winlen
            power(i) = sqrt(mean(Filtchan(i:i+ceil(window*Fs/1000)).^2));
        end
        %detect peaks above detection level (5sd)
        %% remove outliers get better mean value 
        X = power;
        Tempmean(1) = mean(X);
        Tempmean(2) = std(X);
        Xoutlier = X-ones(length(X),1)*Tempmean(1)'...
                            -2.58.*ones(length(X),1)*Tempmean(2)'>0; 
        index = Xoutlier == 1;
        X1 = X;
        X1(index)=[];
        Xfrac = prctile(X1,90);
        X1 = X1.*(X1<Xfrac);
        X1 (X1 == 0) = [];
        Xfrac = prctile(X1,10);
        X1 = X1.*(X1>Xfrac);
        X1 (X1 == 0) = [];
        meanpower = mean(X1,1);% mean
        stdpower = std(X1,0,1);%

        clear X X1 Xfrac index Tempmean Xoutlier 
        %%
        cross = power > meanpower+detectionlevel*stdpower;
        crossbd = [];
        crossbd(:,1) = find(cross-circshift(cross,1)==1);
        crossbd(:,2) = find(cross-circshift(cross,1)==-1);
        clear cross
        gamma = zeros(size(crossbd,1),3);%begin,max,end,frequency
        for i=1:size(crossbd,1)
            gamma(i,2) = find(power(crossbd(i,1):crossbd(i,2))== max(power(crossbd(i,1):crossbd(i,2))))+crossbd(i,1);
            if gamma(i,2) >= 5*Fs % discard the beginning 5 seconds
            %detect beginning boundary
            gamma(i,1) = find(power(1:gamma(i,2))<(meanpower+boundarylevel*stdpower),1,'last');
            gamma(i,3) = find(power(gamma(i,2):end)<(meanpower+boundarylevel*stdpower),1,'first')+gamma(i,2); 
            %detect end boundary      
            end
        end
        gamma = gamma (gamma(:,1)~=0,:);
        gamma = gamma + ceil((window*Fs/1000)/2); %shift values with half window width

        % delete if duplicate
        for i = size(gamma,1):-1:2
            if (gamma(i,1)== gamma(i-1,1))
                gamma(i,:)=[];
            end
        end
        fprintf(['epoch' num2str(p) '-' gammaname '-ch' num2str(ch) ' Events detected: ' num2str(size(gamma,1)) '\n'])
        clear crossbd 
        %% delete wrong duration gamma>5 sample point
        gamma = gamma(gamma(:,3)-gamma(:,1)>5,:);
        %% detect frequency and remove 50 Hz, pick up at least three cycles
        nEvents = size(gamma,1);
        if ~(nEvents == 0)
            for i = 1:nEvents
                A1 = gamma(i,[1,3]);
                [peaks,locs] = findpeaks(Filtchan(A1(1):A1(2)));
                locs1 = [0;locs(1:end-1)];
                med = locs-locs1;
                fre = 1/(median(med(2:end))/Fs);%Hz
                gamma(i,4)= fre;
                gamma(i,5)= length(peaks);
                clear locs locs1 med fre A1
            end
            gamma = gamma(gamma(:,4)~=50,:);% remove 50 Hz
            fprintf([' after remove 50 Hz: ' num2str(size(gamma,1)) '\n'])
            gamma = gamma(gamma(:,5) >= 3,:);% at least three cycles
            fprintf([' larger than 3 three cycles: ' num2str(size(gamma,1)) '\n'])
        end
        %% remove frequency abroad the bundary of the band        
        if ~(nEvents == 0)
           gamma = gamma(gamma(:,4) >= cutoffhighpass & gamma(:,4)<=cutofflowpass,:);
           fprintf([' after remove outlier Hz: ' num2str(size(gamma,1)) '\n'])
        end            
        gammaevents  = gamma;         
               
        GammaEvents{p}{ch} = gammaevents;%epoch,ch*gamma
        clear stdpower meanpower i Filtchan gamma 
        clear x y t gammaevents
    end    
end
        save([savename '_GammaEvents.mat'],'GammaEvents','detchan','-v7.3')

end