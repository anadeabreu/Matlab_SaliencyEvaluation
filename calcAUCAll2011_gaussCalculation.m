
% FYI
% New is just global 
% NewLocal is just Local
% NewWithLocal is Local + Global


clear all
% close all
clc
myTic = tic;
% allGauss=40;%[0.01 .5 1 2 3 5 7.5 10 12.5 15 17.5 20];

allGauss=0.01:0.01:0.13;

% allGauss=0.08:0.01:0.13;
% allGauss=0.01:0.01:0.07; % BORJIIIIIIIIIIIIIIII

datasetNo = 4; % 1 for Bruce ; 2 CRCNS ; 3 for Judd , 4 for Kootstra and 5 for NUSEF

% addpath(genpath('/lab/raid/models/'))

switch datasetNo

% mapSize = size(setTotalFix);

% curSMap = imresize(curSMap, mapSize, 'bilinear');
% kSize = mapSize(2)*sigma;
% curSMap = imfilter(curSMap, fspecial('gaussian', round([kSize, kSize]*4), kSize));

    case 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bruce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bruceScores = cell(3,length(allGauss));
for gaussSize=1:length(allGauss)



%%% Gaussian blob
dims = [100 100]; sigma = 10; % gaussian standard deviation in pixels <<<<<<<<<<<<<
P = [round(dims(1)/2) round(dims(2)/2)];
[Xm Ym] = meshgrid(-P(2):P(2), -P(1):P(1)); s = sigma ;
gauss = exp(-((( Xm.^2)+( Ym.^2)) ./ (2* s^2)));


% load SWalther once here
dims = [511 681];
% disk = strel('disk',50);

load('/lab/tmpir/u/dicky/extras/BruceDataFixations/origfixdata.mat')

path = '/lab/raid/models/OutputSaliency/';
imgPath = '/lab/tmpir/u/dicky/extras/BruceData/';
models = dir(path);

count = 1;
for modelIndex=3:size(models,1)%+2
    
    if ~strcmp(models(modelIndex).name,'SNewLocal')  &&...
         ~strcmp(models(modelIndex).name,'SNewLG')  &&...
         ~strcmp(models(modelIndex).name,'SNewGlobal')                   
%         ~strcmp(models(modelIndex).name,'SHouNIPS') && ...
%         ~strcmp(models(modelIndex).name,'SHouCVPR') && ...
%         ~strcmp(models(modelIndex).name,'SAIM') && ...
%         ~strcmp(models(modelIndex).name,'SSUN') && ...
%         ~strcmp(models(modelIndex).name,'SITTI') && ...
%         ~strcmp(models(modelIndex).name,'SAWS') && ...
%         ~strcmp(models(modelIndex).name,'SSurprise') && ...
%         ~strcmp(models(modelIndex).name,'SGBVS') && ...
%         ~strcmp(models(modelIndex).name,'SSDDR') && ...        
%         ~strcmp(models(modelIndex).name,'SJudd') && ...        
%         ~strcmp(models(modelIndex).name,'SPQFT') 
        continue;
    end
%     
    if modelIndex == size(models,1)+1
        myModelName = 'Gauss';
    elseif modelIndex == size(models,1)+2
        myModelName = 'Human';
        load('/lab/raid/models/AIM/eyetrackingdata/fixdens/allWhite.mat');
    else
        myModelName = models(modelIndex).name;
    end
    
    if strcmp(myModelName,'SAchanta'); continue; end
    if strcmp(myModelName,'SEDS'); continue; end
    if strcmp(myModelName,'SITTI98'); continue; end
%     if strcmp(myModelName,'SHamed'); continue; end
    
    fprintf('%6ds | Calculating auc for model: %s with gauss %1.2f\n',round(toc(myTic)), myModelName,allGauss(gaussSize));
    
    path1 = [path myModelName];
    if exist(path1,'dir')
        cd(path1);
    end
    %     list = dir(path1);
    coefficients = nan(1,120);
%     row = 1; hum=[]; ran=[];
    humCell=cell(120,1);
    ranCell=cell(120,1);
    
    for i=1:120
        %                 load human map
        %         h = find(list(i).name=='.');
        imageName = num2str(i);
        
        img = imread([imgPath imageName '.jpg']);
        %imwrite(img, ['/lab/raid/models/OutputResults/ImagesBruce/' num2str(i) '/' imageName '.jpg'])
        
        clear salMap;
        
        switch myModelName
            case 'SAvraham'
                imageName = ['Esaliencymap' imageName];
            case 'SLeMeur'
                imageName = [imageName 'DensityMap'];
            case 'SJiaLi'
                imageName = ['Jack' num2str(i,'%.6d') '_Sal'];
            case {'SITTI','SITTI98','SSurprise','SEntropy','SVariance'}
                imageName = [imageName '.jpg-VCO000000'];
                
        end
        
        
        switch myModelName
            case {'SLeMeur','SJiaLi','SYan','SGert', 'SNewLocal' ,'SNewLG', 'SNewGlobal' }
                if exist([imageName '.jpg'],'file')
                    salMap = imread([imageName '.jpg']);
                    salMap = salMap(:,:,1);
                elseif exist([imageName '.png'],'file')
                    salMap = imread([imageName '.png']);
                    salMap = salMap(:,:,1);
                end
                
                
            case 'SAWS'
                salMap = load([imageName '.mat']);
                salMap = salMap.SaliencyMap;
                
            case {'SRarityLocal','SRarityGlobal'}
                salMap = load([imageName '.mat']);
                salMap = salMap.att;
                
            case 'SHouCVPR'
                salMap = load([imageName '.mat']);
                salMap = salMap.saliencyMap;
                
            case 'SAIM'
                salMap = load([imageName '.mat']);
                salMap = salMap.info;

%             case 'SSVM'
%                  salMap = load([imageName '.mat']);
%                  salMap = salMap.predictions;
                
            case {'SSVM'}
                salMap = load([imageName '.mat']);
                salMap = salMap.predictions;
            case {'SHamed'}
                salMap = load([imageName '.mat']);
                salMap = salMap.saliency;

            case {'SHouNIPS','SITTI','SITTI98','SVocus','SAvraham','SSurprise','SEntropy','SVariance'}
                if exist([imageName '.jpg'],'file')
                    salMap = double(imread([imageName '.jpg']));
                elseif exist([imageName '.png'],'file')
                    salMap = imread([imageName '.png']);
                end
                salMap = salMap(:,:,1);
                %             case 'SEDS'
                %                 salMap = imread(list(i).name);
                %                 salMap = 255 - salMap;
                %
            case 'SGBVS'
                salMap = load([imageName '.mat']);
                salMap = salMap.master_map;
                
            case 'SJudd'
                salMap = load([imageName '.mat']);
                salMap = salMap.k;
                
            case 'SPQFT'
                salMap = load([imageName '.mat']);
                salMap = salMap.sM;
                
            case 'SSDDR'
                salMap = load([imageName '.mat']);
                salMap = salMap.smap;
                
            case 'SSUN'
                salMap = load([imageName '.mat']);
                salMap = salMap.simg2;
                
            case 'STorallba'
                salMap = load([imageName '.mat']);
                salMap = salMap.map;
                
            case 'SYinLi'
                salMap = load([imageName '.mat']);
                salMap = salMap.mySMap;
                
            case 'SMarat'
                salMap = load(['s_' imageName '.mat']);
                salMap = salMap.S;
                
            case 'SBian'
                salMap = load([imageName '.mat']);
                salMap = salMap.SMAP;
                
            case 'SWalther'
                salMap = load([imageName '.mat']);
                salMap = salMap.sM.data;
                
            case 'Gauss'
                salMap = gauss;
                
            case 'Human'
                %salMap = load([imageName '.mat']);
%                 myIdx = 1:120;
%                 myIdx = myIdx(myIdx ~= i);
%                 salMap = zeros(size(white{i}));
%                 for zz=myIdx
%                     salMap = salMap + white{zz};
%                 end
%                 salMap = conv2(salMap,gauss,'same');
                 salMap = conv2(white{i},gauss,'same');
                
                
        end
        
        %%%% Human
        if strcmp(myModelName,'Human')
            idxNo = i;
            %humanData.eyeMap;
            tmpAC = nan(size(allWhite,1),1);
            for subjNo = 1:size(allWhite,1)
                humanMap = allWhite{subjNo,idxNo};
                
                salMap = zeros(size(humanMap));
                for xyz = 1:size(allWhite,1)
                    if xyz == subjNo; continue; end
                    salMap = salMap + allWhite{xyz,idxNo};
                end
                salMap = conv2(salMap,gauss,'same');
                
                %%% AUC part %%%%%%%%%%
                % Eye fix vs shuffled eye fix of other img
                [X Y] = find(humanMap > 0);
                if isempty(X); continue;end
                idx = 1:size(allWhite,1); idx = idx(idx ~= i);
                randIdx = randperm(length(idx));
                idx = idx(randIdx(1:10)); % only take random 10 pics instead of all 120
                otherHumanMap = zeros(size(humanMap));
                for zidx=idx
                    %filed=strcat('/lab/tmpir/u/dicky/extras/Judd/ALLFIXATIONMAPS/',itsAllFixFiles(zidx+2).name);
                    %fixItem=open(filed);
                    otherHumanMap = otherHumanMap + white{zidx};
                end
                [XRest YRest] = find(otherHumanMap ~= 0);
                
                % Randomized pick
                localHum = nan(length(X),1);
                localRan = nan(length(X),100);
                for k=1:length(X)
                    localHum(k,1) = salMap(X(k),Y(k)); % be aware
                    for kk=1:100
                        r = randi([1 length(XRest)],1);
                        localRan(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                    end
                    %row = row + 1;
                end
                % scoring
                [acz,R] = getauc(localHum, localRan, '');
                tmpAC(subjNo) = mean(acz);
            end
            ac = tmpAC(~isnan(tmpAC));
            coefficients(i) = mean(ac);
        else
        
        
        %%%% Non-human
        if exist('salMap') % for checking SGert
            humanMap = white {i};
%             humanMap = conv2(humanMap,gauss,'same');
            
%             humanMap = (humanMap - min(min(humanMap)))/ (max(max(humanMap)) - min(min(humanMap)));
            
            salMap = double(imresize(salMap,dims,'nearest'));
            
            %%% Gaussing
            mapSize = size(salMap);
            kSize = mapSize(2)*allGauss(gaussSize);
            salMap = imfilter(salMap, fspecial('gaussian', round([kSize, kSize]*4), kSize));
%             salMap = conv2(salMap,gauss,'same');


            if (max(max(salMap)) - min(min(salMap))) == 0;
%                 r = NaN;
                continue
            else
            
            salMap = salMap - min(min(salMap));
            salMap = salMap / max(max(salMap));
            

            %%% AUC part %%%%%%%%%%
            [X Y] = find(humanMap > 0);
            idx = 1:120; idx = idx(idx ~= i);
            otherHumanMap = zeros(size(humanMap));
            for zidx=idx
                otherHumanMap = otherHumanMap + white{zidx};
            end
            [XRest YRest] = find(otherHumanMap ~= 0);
%             X = sacData(find(sacData(:,1)==imgIndex),2);
%             Y = sacData(find(sacData(:,1)==imgIndex),3);
%             
%             XRest = sacData(find(sacData(:,1)~=imgIndex),2);
%             YRest = sacData(find(sacData(:,1)~=imgIndex),3);
            localHum = nan(length(X),1);
            localRan = nan(length(X),100);
            for k=1:length(X)
%                 hum(row,1) = salMap(X(k),Y(k)); % be aware
                localHum(k,1) = salMap(X(k),Y(k)); 
                for kk=1:100
                    r = randi([1 length(XRest)],1);
%                     ran(row,kk) = salMap(XRest(r),YRest(r)); % be aware
                    localRan(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                end
%                 row = row + 1;
            end
            humCell{i} = localHum;
            ranCell{i} = localRan;
            addpath('/lab/raid/models/OutputResults')
            [ac,R] = getauc(localHum, localRan, '');
            
            %%%%%%%%%%%%%%%%%%%%%%%
%             r = corr2(salMap, humanMap);
            coefficients(i) = mean(ac);
            end
            %             [ttr ttp] = corrcoef(salMap,humanMap);
%         else
%             r  = NaN;
        end
        end
%         coefficients = [coefficients r];
        fprintf('.')
    end
    
    
%     [ac,R] = getauc(hum, ran, '');
    %     correlations(count).mean = mean(coefficients);
    %     correlations(count).std = std(coefficients);
    %     correlations(count).coef = coefficients;
    %     correlations(count).name = myModelName;
    AUCValues(count).mean = mean(coefficients);
    AUCValues(count).std = std(coefficients);
%     AUCValues(count).ttp = ttp;
%     AUCValues(count).R = R;
    AUCValues(count).hum=humCell;
    AUCValues(count).ran=ranCell;
%     AUCValues(count).ac=ac;
    AUCValues(count).name = myModelName;
    AUCValues(count).coef = coefficients;
    count = count + 1;
    clear hum ran ac R;
    fprintf('\n')
end

%%% AUC bar chart
myNum = ones(1,length(AUCValues)) * 120;

for ii=1:length(AUCValues)
    tmp = sum(~isnan(AUCValues(ii).coef));
    %if tmp ~= myNum(ii); fprintf('Model %s has only %d out of %d\n', AUCValues(ii).name, tmp, myNum(ii)); end
    myNum(ii) = tmp;
    AUCValues(ii).mean = mean(AUCValues(ii).coef(~isnan(AUCValues(ii).coef)));
    AUCValues(ii).std = std(AUCValues(ii).coef(~isnan(AUCValues(ii).coef)));
end


% plot
xx = vertcat(AUCValues.mean);
yy = vertcat(AUCValues.std) ./ sqrt(myNum');
nn = {AUCValues.name};



cd('/lab/raid/models/OutputResults');

bruceScores{1,gaussSize} = xx;
bruceScores{2,gaussSize} = yy;
bruceScores{3,gaussSize} = nn; % Legend

end

cd('/lab/raid/models/OutputResults');
save('./NewBottom-UpModel/Nov13/BruceLG.mat','-v7.3','bruceScores')

% testgsegsegs
% doesnt work
%errorbar_region(allGauss,horzcat(bruceScores{1,:})',horzcat(bruceScores{2,:})');
% working ver
%plot(allGauss,horzcat(bruceScores{1,:}))


case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CRCNS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


crcnsScores = cell(3,length(allGauss));
for gaussSize=1:length(allGauss)



% myNumRandomSample = 10;
myNumRandomSample = 100;

%%% Gaussian blob
dims = [50 50]; sigma = allGauss(gaussSize); % gaussian standard deviation in pixels <<<<<<<<<<<<<
P = [round(dims(1)/2) round(dims(2)/2)];
[Xm Ym] = meshgrid(-P(2):P(2), -P(1):P(1)); s = sigma ;
gauss = exp(-((( Xm.^2)+( Ym.^2)) ./ (2* s^2)));

%%% all the CRCNS dataset movie names
fldrlist = { ...
    'beverly01'
    'beverly03'
    'beverly05'
    'beverly06'
    'beverly07'
    'beverly08'
    'gamecube02'
    'gamecube04'
    'gamecube05'
    'gamecube06'
    'gamecube13'
    'gamecube16'
    'gamecube17'
    'gamecube18'
    'gamecube23'
    'monica03'
    'monica04'
    'monica05'
    'monica06'
    'saccadetest'
    'standard01'
    'standard02'
    'standard03'
    'standard04'
    'standard05'
    'standard06'
    'standard07'
    'tv-action01'
    'tv-ads01'
    'tv-ads02'
    'tv-ads03'
    'tv-ads04'
    'tv-announce01'
    'tv-music01'
    'tv-news01'
    'tv-news02'
    'tv-news03'
    'tv-news04'
    'tv-news05'
    'tv-news06'
    'tv-news09'
    'tv-sports01'
    'tv-sports02'
    'tv-sports03'
    'tv-sports04'
    'tv-sports05'
    'tv-talk01'
    'tv-talk03'
    'tv-talk04'
    'tv-talk05'
    };

myEyePath = '/lab/raid/models/OutputResults/CRCNSdata/';

%%% For pfmread
addpath('/lab/tmpir/u/dicky/newsal/saliency/matlab');
addpath('/lab/raid/models/OutputResults')


%path2 = '/lab/raid/models/NoBiasData/Judd/ALLFIXATIONMAPS/';
path = '/lab/raid/models/OutputSaliencyCRCNS/';

% itsAllFixFiles = dir('/lab/tmpir/u/dicky/extras/Judd/ALLFIXATIONMAPS/');
% itsAllStiFiles = dir('/lab/tmpir/u/dicky/extras/Judd/ALLSTIMULI/');
% itsNumFixFiles = length(itsAllFixFiles);

models = dir(path);

count = 1;
for modelIndex=3:size(models,1)
%for modelIndex=size(models,1)+1:size(models,1)+2
    
    if modelIndex == size(models,1)+1
        myModelName = 'Gauss';
    elseif modelIndex == size(models,1)+2
        myModelName = 'Human';
    else
        myModelName = models(modelIndex).name;
    end
    
    if strcmp(myModelName,'SAchanta'); continue; end
    if strcmp(myModelName,'SEDS'); continue; end
    if strcmp(myModelName,'SMarat2'); continue; end
    if strcmp(myModelName,'SVarianceCIO'); continue; end
    if strcmp(myModelName,'SEntropyCIO'); continue; end
    if strcmp(myModelName,'SEntropyIgnore'); continue; end
    %     if strcmp(myModelName,'SJiaLi'); continue; end
    
    fprintf('%6ds | Calculating cc for model: %s \n', round(toc(myTic)) ,myModelName);
    
    
    path1 = [path myModelName];
    if strcmp(myModelName,'SMarat') % Exception for SMarat
        path1 = '/lab/raid/models/Sophie/saliency_video';
    end
    if exist(path1,'dir')
        cd(path1);
        %list = dir(path1);
    end
    
    tmpz = dir('./');
    if length(tmpz) ~= length(fldrlist)+2
        fprintf('%s folder list is not equal (%d vs %d)\n',myModelName,length(tmpz)-2,length(fldrlist))
    end
    
    ceofCell = cell(1,length(fldrlist));
    for i=1:length(fldrlist)
        
        %%% Load eye data
        load([myEyePath fldrlist{i} '.mat'])
        
        %%% Check frame length
        frameMax = 999999;
        for subjNo=1:length(eyePos)
            if ~isempty(eyePos{subjNo})
                frameMax = min(frameMax,size(eyePos{subjNo},1));
            end
        end
        
        %%% List the file names & flags
        type='png'; list = dir([fldrlist{i} '.mpg/*png']);
        if strcmp(myModelName,'SMarat'); type='mat'; list = dir([fldrlist{i} '/*mat']); end
        if strcmp(myModelName,'SVOCUS'); type='png'; list = dir(['out' fldrlist{i} '.tgz/*png']); end
        if isempty(list); type='pfm'; list = dir([fldrlist{i} '.mpg/*pfm']); end
        if isempty(list); type='mat'; list = dir([fldrlist{i} '.mpg/*mat']); end
        if strcmp(myModelName,'Gauss') || strcmp(myModelName, 'Human')
        else
            if isempty(list); error('type unknown'); end
        end
        
        
        %%% The index of non-missing subject eyemap
        eyeMapIdx = find(~cellfun('isempty',eyeMap));
        
        %%% Some pretty function
        percNum = 4;
        percNumList = floor(frameMax/percNum:frameMax/percNum:frameMax);
        
        %%% The framemax checking is after nan making due to ttest problem
        coefFrame = nan(frameMax,1);
        
        %%% The file might not be equal to the frame
        if strcmp(myModelName,'Gauss') || strcmp(myModelName, 'Human')
        else
            frameMax = min(frameMax, length(list));
        end
        
        for frameNo = 1:frameMax
            %if max(percNumList == frameNo)
                fprintf('%6ds | Model %s Folder %.2d/%.2d Frame %.4d/%.4d (%3d%%)\n', ...
                    round(toc(myTic)), myModelName, i, length(fldrlist), frameNo, frameMax,round(frameNo*100/frameMax));
            %end
            
            
            if strcmp(myModelName,'Human')
                mepScore = [];
                for leSubj=eyeMapIdx
                    humanMap = zeros(480,640);
                    for zi=eyeMapIdx
                        if zi==leSubj; continue; end
                        humanMap = humanMap + eyeMap{zi}{frameNo};
                    end
                    salMap = conv2(humanMap,gauss,'same');
                    salMap = salMap - min(min(salMap));
                    salMap = salMap / max(max(salMap));
                    
                    humanMap = eyeMap{leSubj}{frameNo};
                    
                    %%% Calculate AUC
                    [X Y] = find(humanMap > 0);
                    idx = 1:frameMax; idx = idx(idx ~= frameNo);
                    randIdx = randperm(length(idx));
                    idx = idx(randIdx(1:10)); % only take random 10 pics instead of all
                    otherHumanMap = zeros(size(humanMap));
                    for zidx=idx
                        for zi=eyeMapIdx
                            otherHumanMap = otherHumanMap + eyeMap{zi}{zidx};
                        end
                    end
                    [XRest YRest] = find(otherHumanMap ~= 0);
                    localHum = nan(length(X),1);
                    localRan = nan(length(X),myNumRandomSample);
                    for k=1:length(X)
                        localHum(k,1) = salMap(X(k),Y(k));
                        for kk=1:myNumRandomSample
                            r = randi([1 length(XRest)],1);
                            localRan(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                        end
                    end
                    [ac,R] = getauc(localHum, localRan, '');
                    
                    mepScore = [mepScore ac];
                end
                coefFrame(frameNo) = mean(mepScore);
            else
                
                
                
                
                %%% Load salmap
                switch type
                    case 'png'
                        if strcmp(myModelName,'SVOCUS')
                            salMap = imread(['out' fldrlist{i} '.tgz/' list(frameNo).name]);
                        else
                            salMap = imread([fldrlist{i} '.mpg/' list(frameNo).name]);
                        end
                    case 'pfm'
                        salMap = pfmread([fldrlist{i} '.mpg/' list(frameNo).name]);
                    case 'mat'
                        if strcmp(myModelName,'SMarat')
                            tmp = load([fldrlist{i} '/' list(frameNo).name]);
                            salMap = tmp.fus_sd; % hey, this is in 640x480 meanwhile others are in 30x40
                            salMap = imresize(salMap,[30 40]);
                        else
                            tmp = load([fldrlist{i} '.mpg/' list(frameNo).name]);
                            myFieldNames = fieldnames(tmp);
                            if length(myFieldNames) ~= 1; error('field of salmap isnt one'); end
                            salMap = tmp.(myFieldNames{1});
                        end
                end
                if strcmp(myModelName,'Gauss')
                    salMap = gauss;
                end
                
                %%% Load eyemap for this frame
                humanMap = zeros(480,640);
                for zi=eyeMapIdx
                    humanMap = humanMap + eyeMap{zi}{frameNo};
                end
                %             humanMap = conv2(humanMap,gauss,'same');
                %             humanMap = humanMap  - min(min(humanMap));
                %             humanMap = humanMap  / max(max(humanMap));
                
                %%% Resize and normalize salmap
                salMap = double(imresize(salMap,size(humanMap),'nearest'));
                salMap = conv2(salMap,gauss,'same');
                salMap = salMap - min(min(salMap));
                salMap = salMap / max(max(salMap));
                
                %%% Showing sample
                %figure;
                %subplot(2,1,1);imshow(humanMap)
                %subplot(2,1,2);imshow(salMap)
                
                %%% Calculate AUC
                %             [X Y] = find(humanMap ~= 0);
                %             NSSVector = zeros(1,size(X,1));
                %             for p=1:size(X,1)
                %                 NSSVector(p) = salMap(X(p),Y(p));
                %             end
                %             coefFrame(frameNo) = mean(NSSVector);
                
                
                [X Y] = find(humanMap > 0);
                idx = 1:frameMax; idx = idx(idx ~= frameNo);
                randIdx = randperm(length(idx));
                idx = idx(randIdx(1:10)); % only take random 10 pics instead of all
                otherHumanMap = zeros(size(humanMap));
                for zidx=idx
                    for zi=eyeMapIdx
                        otherHumanMap = otherHumanMap + eyeMap{zi}{zidx};
                    end
                end
                [XRest YRest] = find(otherHumanMap ~= 0);
                localHum = nan(length(X),1);
                localRan = nan(length(X),myNumRandomSample);
                for k=1:length(X)
                    localHum(k,1) = salMap(X(k),Y(k));
                    for kk=1:myNumRandomSample
                        r = randi([1 length(XRest)],1);
                        localRan(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                    end
                end
                [ac,R] = getauc(localHum, localRan, '');
                coefFrame(frameNo) = mean(ac);
            end
            
            
        end
        ceofCell{i} = coefFrame;
    end
    
    %%% Put in a container
    coefficients = vertcat(ceofCell{:});
    correlations(count).mean = mean(coefficients(~isnan(coefficients)));
    correlations(count).std = std(coefficients(~isnan(coefficients)));
    correlations(count).coef = coefficients;
    correlations(count).name = myModelName;
    count = count + 1;
    
    %     fprintf('\n')
    
end
% error('balh')


myNum = nan(1,length(correlations));

for ii=1:length(correlations)
    tmp = sum(~isnan(correlations(ii).coef));
    %     if tmp ~= myNum(ii); fprintf('Model %s has only %d out of %d\n', correlations(ii).name, tmp, myNum(ii)); end
    myNum(ii) = tmp;
end

% plot
xx = vertcat(correlations.mean);
yy = vertcat(correlations.std) ./ sqrt(myNum');
nn = {correlations.name};

crcnsScores{1,gaussSize} = xx;
crcnsScores{2,gaussSize} = yy;
crcnsScores{3,gaussSize} = nn;

end




    case 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Judd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

juddScores = cell(3,length(allGauss));
for gaussSize=1:length(allGauss)
    
    %%% Gaussian blob
dims = [50 50]; sigma = 10; % gaussian standard deviation in pixels <<<<<<<<<<<<<
P = [round(dims(1)/2) round(dims(2)/2)];
[Xm Ym] = meshgrid(-P(2):P(2), -P(1):P(1)); s = sigma ;
gauss = exp(-((( Xm.^2)+( Ym.^2)) ./ (2* s^2)));


% dims = [768 1024];
% disk = strel('disk',50);

path2 = '/lab/raid/models/NoBiasData/Judd/ALLFIXATIONMAPS/';
path = '/lab/raid/models/OutputSaliencyJudd/';

itsAllFixFiles = dir('/lab/tmpir/u/dicky/extras/Judd/ALLFIXATIONMAPS/');
itsAllStiFiles = dir('/lab/tmpir/u/dicky/extras/Judd/ALLSTIMULI/');
itsNumFixFiles = length(itsAllFixFiles);

humanData = load('/lab/raid/models/NoBiasData/Judd/eyeMap.mat');
addpath('/lab/raid/models/OutputResults/')
models = dir(path);

count = 1;
for modelIndex=3:size(models,1)%+2

    if ~strcmp(models(modelIndex).name,'SNewLocal')  && ...
        ~strcmp(models(modelIndex).name,'SNewLG') && ...                   
        ~strcmp(models(modelIndex).name,'SNewGlobal')                   
%         ~strcmp(models(modelIndex).name,'SHouNIPS') && ...
%         ~strcmp(models(modelIndex).name,'SHouCVPR') && ...
%         ~strcmp(models(modelIndex).name,'SAIM') && ...
%         ~strcmp(models(modelIndex).name,'SSUN') && ...
%         ~strcmp(models(modelIndex).name,'SITTI') && ...
%         ~strcmp(models(modelIndex).name,'SAWS') && ...
%         ~strcmp(models(modelIndex).name,'SSurprise') && ...
%         ~strcmp(models(modelIndex).name,'SGBVS') && ...
%         ~strcmp(models(modelIndex).name,'SSDDR') && ...        
%         ~strcmp(models(modelIndex).name,'SJudd') && ...        
%         ~strcmp(models(modelIndex).name,'SPQFT') 
        continue;
    end
    
%     if ~strcmp(models(modelIndex).name,'SNewLocal') && ...
%         ~strcmp(models(modelIndex).name,'SNewWithLocal') 
%         continue;
%     end    
    

    if modelIndex == size(models,1)+1
        myModelName = 'Gauss';
    elseif modelIndex == size(models,1)+2
        myModelName = 'Human';
    else
        myModelName = models(modelIndex).name;
    end
    
    if strcmp(myModelName,'SAchanta'); continue; end
    if strcmp(myModelName,'SEDS'); continue; end
    %if strcmp(myModelName,'SHamed'); continue; end
    if strcmp(myModelName,'SJudd2'); continue; end
%     if strcmp(myModelName,'SJiaLi'); continue; end    
    
    fprintf('calculating AUC for model: %s with gauss %1.2f\n', myModelName,allGauss(gaussSize));
    
    
    path1 = [path myModelName];
    if exist(path1,'dir')
    cd(path1);
    list = dir(path1);
    end
    
    coefficients = nan(1,size(itsAllFixFiles,1)-2);
    row = 1; hum=[]; ran=[];
    for i=3:size(itsAllFixFiles,1)
        
        if strcmp(myModelName, 'SVocus')
            if i==3; continue; end
            i = i-1; %%% Vocus got wrong order images
        end
        
        % load human map
        h = find(itsAllFixFiles(i).name=='.');
        imageName = itsAllFixFiles(i).name(1:h-1);
        
        
        switch myModelName
            case 'SAvraham'
                imageName = ['Esaliencymap' imageName];
            case 'SLeMeur'
                imageName = [imageName 'DensityMap'];
%             case 'SJiaLi'
%                 imageName = ['Jack' num2str(i,'%.6d') '_Sal'];
            case {'SITTI','SITTI98','SSurprise','SEntropy','SVariance'}
                imageName = [imageName '.jpeg-VCO000000'];
                
        end
        
        
        switch myModelName
            case {'SLeMeur','SJiaLi','SYan','SGert', 'SNewLocal','SNewLG','SNewGlobal'}
                if exist([imageName '.jpg'],'file')
                    salMap = imread([imageName '.jpg']);
                elseif exist([imageName '.png'],'file')
                    salMap = imread([imageName '.png']);
                elseif exist([imageName '.jpeg'],'file')
                    salMap = imread([imageName '.jpeg']);
                end
                if exist('salMap','var')
                    salMap = salMap(:,:,1);
                end
                
                
            case 'SAWS'
                salMap = load([imageName '..mat']);
                salMap = salMap.SaliencyMap;
                
            case {'SRarityLocal','SRarityGlobal'}
                salMap = load([imageName '.mat']);
                salMap = salMap.att;
                
            case 'SHouCVPR'
                if exist([imageName '.jpeg.mat'])
                salMap = load([imageName '.jpeg.mat']);
                salMap = salMap.saliencyMap;
                end
                
            case 'SAIM'
                salMap = load([imageName '..mat']);
                salMap = salMap.info;
                
            case {'SITTI','SITTI98','SVocus','SAvraham','SSurprise','SEntropy','SVariance'}
                if exist([imageName '.jpg'],'file')
                    salMap = imread([imageName '.jpg']);
                elseif exist([imageName '.png'],'file')
                    salMap = imread([imageName '.png']);
                elseif exist([imageName '.jpeg'],'file')
                    salMap = imread([imageName '.jpeg']);
                end
                if exist('salMap','var')
                salMap = salMap(:,:,1);
                end
                %             case 'SEDS'
                %                 salMap = imread(list(i).name);
                %                 salMap = 255 - salMap;
                %
            case 'SHouNIPS'
                salMap = load([imageName '..mat']);
                salMap = salMap.mySMap;
            
            case 'SGBVS'
                salMap = load([imageName '.mat']);
                salMap = salMap.master_map;
                
            case 'SJudd'
                salMap = load([imageName '.mat']);
                salMap = salMap.k;
                
            case 'SPQFT'
                salMap = load([imageName '.mat']);
                salMap = salMap.sM;
                
            case 'SSDDR'
                salMap = load([imageName '.mat']);
                salMap = salMap.smap;
                
            case 'SSUN'
                salMap = load([imageName '.mat']);
                salMap = salMap.simg2;
                
            case 'STorallba'
                salMap = load([imageName '.mat']);
                salMap = salMap.map;
                
            case 'SYinLi'
                salMap = load([imageName '..mat']);
                salMap = salMap.mySMap;
                
            case 'SMarat'
                salMap = load(['s_' imageName '.mat']);
                salMap = salMap.S;
                
            case 'SBian'
                salMap = load([imageName '.mat']);
                salMap = salMap.SMAP;
                
            case 'SWalther'
                salMap = load([imageName '.mat']);
                salMap = salMap.sM.data;

            case 'SHamed'
                salMap = load([imageName '.mat']);
                salMap = salMap.finalSal;
            case 'SSVM'
                salMap = load([imageName '.jpeg.mat']);
                salMap = salMap.predictions;
                
            case 'Gauss'
                salMap = gauss;
                
            case 'Human'
                filed=strcat('/lab/tmpir/u/dicky/extras/Judd/ALLFIXATIONMAPS/',itsAllFixFiles(i).name);
                fixItem=open(filed);
                himg=im2double(fixItem.fixationPts);
                salMap = conv2(himg,gauss,'same');
                
        end
        
        if strcmp(myModelName, 'SVocus')
            i = i+1; %%% Vocus got wrong order images; this to revert previous -1
        end
        
        if strcmp(myModelName,'Human')
            idxNo = i-2;
            humanData.eyeMap;
            tmpAC = nan(size(humanData.eyeMap,2),1);
            for subjNo = 1:size(humanData.eyeMap,2)
                humanMap = humanData.eyeMap{idxNo,subjNo};
                
                salMap = zeros(size(humanMap));
                for xyz = 1:size(humanData.eyeMap,2)
                    if xyz == subjNo; continue; end
                    salMap = salMap + humanData.eyeMap{idxNo,xyz};
                end
                salMap = conv2(salMap,gauss,'same');
                
                %%% AUC part %%%%%%%%%%
                % Eye fix vs shuffled eye fix of other img
                [X Y] = find(humanMap > 0);
                idx = 1:size(itsAllFixFiles,1)-2; idx = idx(idx ~= i);
                randIdx = randperm(length(idx));
                idx = idx(randIdx(1:10)); % only take random 10 pics instead of all 1003
                otherHumanMap = zeros(size(humanMap));
                for zidx=idx
                    filed=strcat('/lab/tmpir/u/dicky/extras/Judd/ALLFIXATIONMAPS/',itsAllFixFiles(zidx+2).name);
                    fixItem=open(filed);
                    otherHumanMap = otherHumanMap + imresize(im2double(fixItem.fixationPts),size(otherHumanMap),'nearest');
                end
                [XRest YRest] = find(otherHumanMap ~= 0);
                
                % Randomized pick
                localHum = nan(length(X),1);
                localRan = nan(length(X),100);
                for k=1:length(X)
                    localHum(k,1) = salMap(X(k),Y(k)); % be aware
                    for kk=1:100
                        r = randi([1 length(XRest)],1);
                        localRan(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                    end
                    row = row + 1;
                end
                % scoring
                [acz,R] = getauc(localHum, localRan, '');
                tmpAC(subjNo) = mean(acz);
            end
            ac = tmpAC;
        else
            if ~exist('salMap','var');fprintf('model %s: file %s is missing\n',myModelName,itsAllFixFiles(i).name); continue; end
            
            
            filed=strcat('/lab/tmpir/u/dicky/extras/Judd/ALLFIXATIONMAPS/',itsAllFixFiles(i).name);
            fixItem=open(filed);
            %     himg=im2double(fixItem.saliencyMapBlur);
            humanMap=im2double(fixItem.fixationPts);
            %         humanMap = conv2(humanMap,gauss,'same');
            %         humanMap = humanMap  - min(min(humanMap));
            %         humanMap = humanMap  / max(max(humanMap));
            
            % load model saliency map - AIM
            salMap = double(imresize(salMap,size(humanMap),'nearest'));
            
            %%% Gaussing
            mapSize = size(salMap);
            kSize = mapSize(2)*allGauss(gaussSize);
            salMap = imfilter(salMap, fspecial('gaussian', round([kSize, kSize]*4), kSize));
            
%             salMap = conv2(salMap,gauss,'same');
            salMap = salMap - min(min(salMap));
            salMap = salMap / max(max(salMap));
            %mapMean = mean2(salMap); mapStd = std2(salMap);
            %salMap = (salMap - mapMean) / mapStd; % Normalized map
            
            
            %%% AUC part %%%%%%%%%%
            [X Y] = find(humanMap > 0);
            idx = 1:size(itsAllFixFiles,1)-2; idx = idx(idx ~= i);
            randIdx = randperm(length(idx));
            idx = idx(randIdx(1:10)); % only take random 10 pics instead of all 1003
            otherHumanMap = zeros(size(humanMap));
            for zidx=idx
                filed=strcat('/lab/tmpir/u/dicky/extras/Judd/ALLFIXATIONMAPS/',itsAllFixFiles(zidx+2).name);
                fixItem=open(filed);
                otherHumanMap = otherHumanMap + imresize(im2double(fixItem.fixationPts),size(otherHumanMap),'nearest');
            end
            [XRest YRest] = find(otherHumanMap ~= 0);
            %             X = sacData(find(sacData(:,1)==imgIndex),2);
            %             Y = sacData(find(sacData(:,1)==imgIndex),3);
            %
            %             XRest = sacData(find(sacData(:,1)~=imgIndex),2);
            %             YRest = sacData(find(sacData(:,1)~=imgIndex),3);
            %
            localHum = nan(length(X),1);
            localRan = nan(length(X),100);
            for k=1:length(X)
                %hum(row,1) = salMap(X(k),Y(k)); % be aware
                localHum(k,1) = salMap(X(k),Y(k)); % be aware
                for kk=1:100
                    r = randi([1 length(XRest)],1);
                    %ran(row,kk) = salMap(XRest(r),YRest(r)); % be aware
                    localRan(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                end
                row = row + 1;
            end
            [ac,R] = getauc(localHum, localRan, '');
            
            
            %         [X Y] = find(humanMap ~= 0);
            %         NSSVector = zeros(1,size(X,1));
            %         for p=1:size(X,1)
            %             NSSVector(p) = salMap(X(p),Y(p));
            %         end
            %         r = mean(NSSVector);
            %r = corr2(salMap, humanMap);
            %             [ttr ttp] = corrcoef(salMap,humanMap);
            %         coefficients = [coefficients r];
            
        end
        coefficients(i-2) = mean(ac);
        %fprintf('.')
        if mod(i,100)==0; fprintf('%6ds | %s %.4d/%.4d\n',round(toc(myTic)),myModelName,i,length(list)); end

%         figure; 
%         imshow(salMap,[])
%         pause
        clear salMap;
        
    end
    
% correlations(count).mean = mean(coefficients(~isnan(coefficients)));
% correlations(count).std = std(coefficients(~isnan(coefficients)));
% correlations(count).coef = coefficients;
% correlations(count).name = myModelName;

addpath('/lab/raid/models/OutputResults')
%     [ac,R] = getauc(hum, ran, '');
    %     correlations(count).mean = mean(coefficients);
    %     correlations(count).std = std(coefficients);
    %     correlations(count).coef = coefficients;
    %     correlations(count).name = myModelName;
    AUCValues(count).mean = mean(ac);
    AUCValues(count).std = std(ac);
%     AUCValues(count).ttp = ttp;
%     AUCValues(count).R = R;
    AUCValues(count).ac = coefficients;
%     AUCValues(count).hum=hum;
%     AUCValues(count).ran=ran;
    AUCValues(count).name = myModelName;
    AUCValues(count).coef = coefficients;
%     count = count + 1;
    clear hum ran ac R;
    
count = count + 1;

fprintf('\n')

end
% error('balh')


myNum = ones(1,length(AUCValues)) * (size(itsAllFixFiles,1)-2);


for ii=1:length(AUCValues)
    tmp = sum(~isnan(AUCValues(ii).coef));
    if tmp ~= myNum(ii); fprintf('Model %s has only %d out of %d\n', AUCValues(ii).name, tmp, myNum(ii)); end
    myNum(ii) = tmp;
    AUCValues(ii).mean = mean(AUCValues(ii).coef(~isnan(AUCValues(ii).coef)));
    AUCValues(ii).std = std(AUCValues(ii).coef(~isnan(AUCValues(ii).coef)));
end



xx = vertcat(AUCValues.mean);
yy = vertcat(AUCValues.std) ./ sqrt(myNum');
nn = {AUCValues.name};

cd('/lab/raid/models/OutputResults/NewBottom-UpModel/');
% save(['aucTMPscore_Judd_LG_' num2str(allGauss(gaussSize),'%1.2f') '.mat'],'-v7.3','xx','yy','nn')
    
juddScores{1,gaussSize} = xx;
juddScores{2,gaussSize} = yy;
juddScores{3,gaussSize} = nn;
    
end

cd('/lab/raid/models/OutputResults');
save('./NewBottom-UpModel/Nov13/JuddLG.mat','-v7.3','bruceScores')

% save('./NewBottom-UpModel/aucAllScores_JuddOther.mat','-v7.3','juddScores')
% % % % % 
% % % % % 
% % % % % zzzzgnasdlgnlas

    case 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Kootstra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kootstraScores = cell(3,length(allGauss));
for gaussSize=1:length(allGauss)
    %%% Gaussian blob
dims = [50 50]; sigma = 10; % gaussian standard deviation in pixels <<<<<<<<<<<<<
P = [round(dims(1)/2) round(dims(2)/2)];
[Xm Ym] = meshgrid(-P(2):P(2), -P(1):P(1)); s = sigma ;
gauss = exp(-((( Xm.^2)+( Ym.^2)) ./ (2* s^2)));


% load SWalther once here
dims = [768 1024];
% disk = strel('disk',50);

load('/lab/tmpir/u/dicky/extras/KootstraFixations/eyeTrackData.mat')

path = '/lab/raid/models/OutputSaliencyKootstra/';
models = dir(path);

%%% Kootstra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
types = { ...
    'buildings'
    'nature'
    'animals'
    'flowers'
    'automan'
    };

numbers = { ...
    [0 2 8:21] % buildings
    [1 3 4 6 8 10 11 14:16 18:23 25:35 37:40 42:50] % nature
    [0:11] % animals
    [0:19] % flowers
    [0 2:12] % automan
    };

% making MEP for each subjects and image

% mepmap = zeros(768,1024);
% mepsalmap = zeros(768,1024);
myMEP = cell(1,length(types));
for i=1:length(types)
    myMEPno = cell(1,length(numbers{i}));
    for j=1:length(numbers{i})
        % mep map
        subjs=eyeTrackData.(types{i}).([types{i} '_' num2str(numbers{i}(j),'%.2d')]);
        subjCell = struct2cell(subjs);
        myMEPsubj = cell(1,length(subjCell));
        for k=1:length(subjCell)
            x = subjCell{k}.fixX;
            y = subjCell{k}.fixY;
            if length(x) ~= length(y); error('size mismatch');end
            mepmap = zeros(768,1024);
            for m=1:length(x)
                xx = min(max(round(x(m)),1),1024);
                yy = min(max(round(y(m)),1),768);
                mepmap(yy, xx) = mepmap(yy, xx) + 1;
            end
            myMEPsubj{k} = mepmap;
        end
        myMEPno{j} = myMEPsubj;
        % mepsalmap
        %         if exist([kootstraSalPath types{i} '_' num2str(numbers{i}(j),'%.2d') '.png'],'file')
        %             salmap = imread([kootstraSalPath types{i} '_' num2str(numbers{i}(j),'%.2d') '.png']);
        %             mepsalmap = mepsalmap + double(rgb2gray(salmap));
        %         end
    end
    myMEP{i} = myMEPno;
end
% ----------------------------------


% myMEP{type}{imgnumber}{subj}

count = 1;
correlations = [];

for modelIndex=3:size(models,1)%+2
    
    if ~strcmp(models(modelIndex).name,'SNewLGComb3') 
%         ~strcmp(models(modelIndex).name,'SNewLocal') &&...           
%         ~strcmp(models(modelIndex).name,'SNewGlobal')  
        
        continue;
    
    end
    
    if modelIndex == size(models,1)+1
        myModelName = 'Gauss';
    elseif modelIndex == size(models,1)+2
        myModelName = 'Human'; %%% What???
    else
        myModelName = models(modelIndex).name;
    end
    
    if strcmp(myModelName,'SAchanta'); continue; end
    if strcmp(myModelName,'SEDS'); continue; end
    %     if strcmp(myModelName,'SJiaLi'); continue; end
    
    fprintf('calculating AUC for model: %s with gauss %1.2f\n', myModelName,allGauss(gaussSize));
    
    path1 = [path myModelName];
    if exist(path1,'dir')
        cd(path1);
    end
    coefficients = [];
    row = 1; %hum=[]; ran=[];
%     humCell = cell(1,5);
%     ranCell = cell(1,5);
    
    for catIdx = 1:5
        myCoeffScore = [];
        
%         humC2 = cell(1,length(numbers{catIdx}));
%         ranC2 = cell(1,length(numbers{catIdx}));
        for imgIdx = 1:length(numbers{catIdx})
            imageName = [types{catIdx} '_' num2str(numbers{catIdx}(imgIdx),'%.2d')];%num2str(i);
            clear salMap;
            
            switch myModelName
                case 'SAvraham'
                    imageName = ['Esaliencymap' imageName];
                case 'SLeMeur'
                    imageName = [imageName 'DensityMap'];
                    %case 'SJiaLi'
                    %imageName = ['Jack' num2str(i,'%.6d') '_Sal'];
                case {'SITTI','SITTI98','SSurprise','SEntropy','SVariance'}
                    imageName = [imageName '.png-VCO000000'];
                    
            end
            
            
            switch myModelName
                case {'SLeMeur','SJiaLi','SYan','SGert', 'SNewLocal',  'SNewLGComb3', 'SNewGlobal'}
                    if exist([imageName '.jpg'],'file')
                        salMap = imread([imageName '.jpg']);
                        salMap = salMap(:,:,1);
                    elseif exist([imageName '.png'],'file')
                        salMap = imread([imageName '.png']);
                        salMap = salMap(:,:,1);
                    end
                    
                case 'SAWS'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.SaliencyMap;
                    
                case {'SRarityLocal','SRarityGlobal'}
                    salMap = load([imageName '.mat']);
                    salMap = salMap.att;
                    
                case 'SHouCVPR'
                    salMap = load([imageName '.png.mat']);
                    salMap = salMap.saliencyMap;
                    
                case 'SAIM'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.info;
                    
                case {'SITTI','SITTI98','SVocus','SAvraham','SSurprise','SEntropy','SVariance'}
                    if exist([imageName '.jpg'],'file')
                        salMap = imread([imageName '.jpg']);
                    elseif exist([imageName '.png'],'file')
                        salMap = imread([imageName '.png']);
                    end
                    salMap = salMap(:,:,1);
                    %             case 'SEDS'
                    %                 salMap = imread(list(i).name);
                    %                 salMap = 255 - salMap;
                    %
                case 'SGBVS'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.master_map;
                    
                case 'SJudd'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.k;
                    
                case 'SPQFT'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.sM;
                    
                case 'SSDDR'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.smap;
                    
                case 'SSUN'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.simg2;
                    
                case 'STorallba'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.map;
                    
                case 'SHouNIPS'
                    salMap = load([imageName '.png.mat']);
                    salMap = salMap.mySMap;
                    
                case 'SYinLi'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.mySMap;
                    
                case 'SMarat'
                    salMap = load(['s_' imageName '.mat']);
                    salMap = salMap.S;
                    
                case 'SBian'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.SMAP;
                    
                case 'SWalther'
                    salMap = load([imageName 'mat.mat']);
                    salMap = salMap.sM.data;
                case 'SHamed'
                    salMap = load([imageName '.mat']);
                    salMap = salMap.saliency;
                case 'SSVM'
                    salMap = load([imageName '.png.mat']);
                    if isfield(salMap,'predictions')
                        salMap = salMap.predictions;
                    else
                        salMap = salMap.master_map;
                    end
                    
                case 'Gauss'
                    salMap = gauss;
                    
                    %                 case 'Human'
                    %                     %salMap = load([imageName '.mat']);
                    %                     myIdx = 1:120;
                    %                     myIdx = myIdx(myIdx ~= i);
                    %                     salMap = zeros(size(white{i}));
                    %                     for zz=myIdx
                    %                         salMap = salMap + white{zz};
                    %                     end
                    %                     salMap = conv2(salMap,gauss,'same');
                    
                    
            end
            
            %        figure; imshow(salMap,[])
            if strcmp(myModelName,'Human')
                %humanMap = zeros(768,1024);
                %myMEP{type}{imgnumber}{subj}
%                 humC3 = cell(1,length(myMEP{catIdx}{imgIdx}));
%                 ranC3 = cell(1,length(myMEP{catIdx}{imgIdx}));
                rr = nan(1,length(myMEP{catIdx}{imgIdx}));
                for subjNo = 1:length(myMEP{catIdx}{imgIdx})
                    humanMap = zeros(768,1024);
                    for zzz=1:length(myMEP{catIdx}{imgIdx})
                        if zzz == subjNo; continue; end
                        humanMap = humanMap + myMEP{catIdx}{imgIdx}{zzz};
                    end
                    salMap = conv2(humanMap,gauss,'same');
                    salMap = salMap - min(min(salMap));
                    salMap = salMap / max(max(salMap));
            
                    %humanMap = (humanMap - min(min(humanMap)))/ (max(max(humanMap)) - min(min(humanMap)));
                    %                     salMap = humanMap;
                    %                     mapMean = mean2(salMap); mapStd = std2(salMap);
                    %                     salMap = (salMap - mapMean) / mapStd; % Normalized map
                    
                    realHumanMap = myMEP{catIdx}{imgIdx}{subjNo};
                    
                    %%% AUC
                    [X Y] = find(realHumanMap > 0);
                    
                    %%% We want to use other images
                    idx = 1:length(numbers{catIdx}); idx = idx(idx ~= imgIdx);
                    otherHumanMap = zeros(size(humanMap));
                    for zidx=idx
                        for zzz=1:length(myMEP{catIdx}{zidx})
                            otherHumanMap = otherHumanMap + myMEP{catIdx}{zidx}{zzz};
                        end
                    end
                    [XRest YRest] = find(otherHumanMap ~= 0);
                    
                    hum = zeros(length(X),1);
                    ran = zeros(length(X),100);
                    for k=1:length(X)
                        hum(k,1) = salMap(X(k),Y(k)); % be aware
                        for kk=1:100
                            r = randi([1 length(XRest)],1);
                            ran(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                        end
                        %row = row + 1;
                    end
%                     humC3{subjNo} = hum;
%                     ranC3{subjNo} = ran;
                    
                    [ac,R] = getauc(hum, ran, '');
                    
                    rr(subjNo) = mean(ac);
                    
                    %                     [X Y] = find(realHumanMap ~= 0);
                    %                     NSSVector = zeros(1,size(X,1));
                    %                     for p=1:size(X,1)
                    %                         NSSVector(p) = salMap(X(p),Y(p));
                    %                     end
                    %                     rr{subjNo} = NSSVector;
                    %                     realHumanMap = conv2(realHumanMap,gauss,'same');
                    %                     realHumanMap = (realHumanMap - min(min(realHumanMap)))/ (max(max(realHumanMap)) - min(min(realHumanMap)));
                    
                    %                     rr(subjNo) = corr2(humanMap, realHumanMap);
                end
%                 humC2{imgIdx} = vertcat(humC3{:});
%                 ranC2{imgIdx} = vertcat(ranC3{:});
                %                 r = mean(horzcat(rr{:}));
                r = mean(rr);
                
            else
                if exist('salMap','var') % for checking SGert
                    humanMap = zeros(768,1024);
                    for zzz=1:length(myMEP{catIdx}{imgIdx})
                        humanMap = humanMap + myMEP{catIdx}{imgIdx}{zzz};
                    end
                    
                    %humanMap = conv2(humanMap,gauss,'same');
                    %humanMap = (humanMap - min(min(humanMap)))/ (max(max(humanMap)) - min(min(humanMap)));
                    
                    salMap = double(imresize(salMap,dims,'nearest'));
                    
                    %%% Gaussing
                    mapSize = size(salMap);
                    kSize = mapSize(2)*allGauss(gaussSize);
                    salMap = imfilter(salMap, fspecial('gaussian', round([kSize, kSize]*4), kSize));
%                     salMap = conv2(salMap,gauss,'same');
                    if (max(max(salMap)) - min(min(salMap))) == 0
                        r = NaN;
                    else
                        salMap = salMap - min(min(salMap));
                        salMap = salMap / max(max(salMap));
            
                        
                        %%% For image-difference
                        idx = 1:length(numbers{catIdx}); idx = idx(idx ~= imgIdx);
                        otherHumanMap = zeros(size(humanMap));
                        for zidx=idx
                            for zzz=1:length(myMEP{catIdx}{zidx})
                                otherHumanMap = otherHumanMap + myMEP{catIdx}{zidx}{zzz};
                            end
                        end
                        [X Y] = find(humanMap > 0);
                        [XRest YRest] = find(otherHumanMap > 0);
                        hum = zeros(length(X),1);
                        ran = zeros(length(X),100);
                        for k=1:length(X)
                            hum(k,1) = salMap(X(k),Y(k)); % be aware
                            for kk=1:100
                                r = randi([1 length(XRest)],1);
                                ran(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                            end
                        end
                        [ac,R] = getauc(hum, ran, '');
%                         humC2{imgIdx} = hum;
%                         ranC2{imgIdx} = ran;
                        
                        r = mean(ac);
                    end
                else
                    r  = NaN;
                end
            end
            coefficients = [coefficients r];
            myCoeffScore = [myCoeffScore r];
            %             fprintf('.')
            fprintf('%6ds | %s %.1d/%.1d -> %.2d/%.2d\n',round(toc(myTic)),myModelName,catIdx,length(types),imgIdx,length(numbers{catIdx}));
        end
        %         correlationsCat(count,catIdx).mean = mean(myCoeffScore);
        %         correlationsCat(count,catIdx).std = std(myCoeffScore);
        %         correlationsCat(count,catIdx).coef = myCoeffScore;
        %         correlationsCat(count,catIdx).name = myModelName;
        %         correlationsCat(count,catIdx).catname =types{catIdx};
%         hum = vertcat(humC2{:});
%         ran = vertcat(ranC2{:});
%         [ac,R] = getauc(hum, ran, '');
        AUCValuesCat(count,catIdx).mean = mean(myCoeffScore);
        AUCValuesCat(count,catIdx).std = std(myCoeffScore);
%         AUCValuesCat(count,catIdx).R = R;
%         AUCValuesCat(count,catIdx).hum=hum;
%         AUCValuesCat(count,catIdx).ran=ran;
        AUCValuesCat(count,catIdx).ac=myCoeffScore;
        AUCValuesCat(count,catIdx).name = myModelName;
        AUCValuesCat(count,catIdx).coef = myCoeffScore;
        AUCValuesCat(count,catIdx).catname = types{catIdx};
%         humCell{catIdx} = humC2;
%         ranCell{catIdx} = ranC2;
    end
%     hum = vertcat(humCell{1}{:}, humCell{2}{:}, humCell{3}{:}, humCell{4}{:}, humCell{5}{:});
%     ran = vertcat(ranCell{1}{:}, ranCell{2}{:}, ranCell{3}{:}, ranCell{4}{:}, ranCell{5}{:});
%     [ac,R] = getauc(hum, ran, '');
    AUCValues(count).mean = mean(coefficients);
    AUCValues(count).std = std(coefficients);
%     AUCValues(count).R = R;
%     AUCValues(count).hum=hum;
%     AUCValues(count).ran=ran;
    AUCValues(count).ac=coefficients;
    AUCValues(count).name = myModelName;
    AUCValues(count).coef = coefficients;
    count = count + 1;
    
    %     correlations(count).mean = mean(coefficients);
    %     correlations(count).std = std(coefficients);
    %     correlations(count).coef = coefficients;
    %     correlations(count).name = myModelName;
    %     count = count + 1;
    fprintf('\n')
end
%
% cd('/lab/raid/models/OutputResults');
% save('ccBruceNew2011.mat','correlations')

% error('blah')

myNum = ones(1,length(AUCValues)) * 100;


for zing=1:length(AUCValues)
    ii = ~isnan(AUCValues(zing).ac);
    AUCValues(zing).mean = mean(AUCValues(zing).ac(ii));
    AUCValues(zing).std = std(AUCValues(zing).ac(ii));
    myNum(zing) = sum(ii);
end

%%% ploting %%%%%%%%%%%%%%%%%%%%%%%%%
xx = vertcat(AUCValues.mean);
yy = vertcat(AUCValues.std) ./ sqrt(myNum');
nn = {AUCValues.name};


kootstraScores{1,gaussSize} = xx;
kootstraScores{2,gaussSize} = yy;
kootstraScores{3,gaussSize} = nn;
    
end

cd('/lab/raid/models/OutputResults');
save('./NewBottom-UpModel/Nov13/KootstraLGComb.mat','-v7.3','kootstraScores')
% % % % % % 
% % % % % % hadmgklsgdlsjldgjsl


    case 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% nusef
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nusefScores = cell(3,length(allGauss));
for gaussSize=1:length(allGauss)
    
    %%% Gaussian blob
dims = [50 50]; sigma = 10; % gaussian standard deviation in pixels <<<<<<<<<<<<<
P = [round(dims(1)/2) round(dims(2)/2)];
[Xm Ym] = meshgrid(-P(2):P(2), -P(1):P(1)); s = sigma ;
gauss = exp(-((( Xm.^2)+( Ym.^2)) ./ (2* s^2)));

addpath('/lab/raid/models/OutputResults/')
addpath('/lab/tmpir/u/dicky/nusef/NUSEF_database/code/')

path = '/lab/tmpir/u/dicky/nusef/OutputSaliency/';

itsAllStiFiles = dir('/lab/tmpir/u/dicky/nusef/NUSEF_database/stimuli/*.jpg');
itsNumFixFiles = length(itsAllStiFiles);

models = dir(path);

count = 1;
for modelIndex=3:size(models,1)%+2%2
    
     if       ~strcmp(models(modelIndex).name,'SNewLocal') && ...
              ~strcmp(models(modelIndex).name,'SNewLG') && ...              
              ~strcmp(models(modelIndex).name,'SNewGlobal') 
%             ~strcmp(models(modelIndex).name,'SHouCVPR') && ...
%             ~strcmp(models(modelIndex).name,'SAIM') && ...
%             ~strcmp(models(modelIndex).name,'SSUN') && ...
%             ~strcmp(models(modelIndex).name,'SITTI') && ...
%             ~strcmp(models(modelIndex).name,'SAWS') && ...
%             ~strcmp(models(modelIndex).name,'SSurprise') && ...
%             ~strcmp(models(modelIndex).name,'SGBVS') && ...
%             ~strcmp(models(modelIndex).name,'SSDDR') && ...
%             ~strcmp(models(modelIndex).name,'SJudd') && ...
%             ~strcmp(models(modelIndex).name,'SPQFT')
         continue;
     end
   

%     if ~strcmp(models(modelIndex).name,'SNewLocal') &&...
%             ~strcmp(models(modelIndex).name,'SNewWithLocal') 
%          continue;
%         
%     end

    
    if modelIndex == size(models,1)+1
        myModelName = 'Gauss';
    elseif modelIndex == size(models,1)+2
        myModelName = 'Human';
    else
        myModelName = models(modelIndex).name;
    end
    
    fprintf('calculating AUC for model: %s with gauss %1.2f\n', myModelName, allGauss(gaussSize));
    
    
    path1 = [path myModelName];
    if exist(path1,'dir')
        cd(path1);
        list = dir(path1);
    end
    
    coefficients = nan(1,size(itsAllStiFiles,1));
    row = 1; hum=[]; ran=[];
    for i=1:size(itsAllStiFiles,1)
        
        % The image name
        h = find(itsAllStiFiles(i).name=='.');
        imageName = itsAllStiFiles(i).name(1:h-1);
        
        if 1==2%exist(['/lab/tmpir/u/dicky/nusef/cache/' myModelName '_' imageName '_' num2str(allGauss(gaussSize),'%1.2f') '.mat'],'file')
             load(['/lab/tmpir/u/dicky/nusef/cache/' myModelName '_' imageName '_' num2str(allGauss(gaussSize),'%1.2f') '.mat'])
        else
            
            switch myModelName
                case 'Gauss'
                    salMap = gauss;
                case 'Human'
                    salMap = [];
                case 'SNewLocal'
                    load([imageName '.mat'])
                    salMap = salMapL;
                case 'SNewGlobal'
                    load([imageName '.mat'])
                    salMap = salMapG;
                case 'SNewLG'
                    load([imageName '.mat'])
                    salMap = salMapLG;
                    
                otherwise
                    load([imageName '.mat'])
                    salMap = SM;
                    clear SM
            end
            
            if ~exist('salMap','var');fprintf('model %s: file %s is missing\n',myModelName,itsAllStiFiles(i).name); continue; end
            fixItem = get_NUSEF_fix_data(itsAllStiFiles(i).name);
            
            if strcmp(myModelName,'Human')
                idxSubj = 1:length(fixItem.fix);
                tmpAC = nan(1,length(fixItem.fix));
                for subjNo=idxSubj
                    idxOther = idxSubj(idxSubj~=subjNo);
                    humanMap = zeros(size(fixItem.map));
                    salMap = zeros(size(fixItem.map));
                    
                    humanFix = fixItem.fix{subjNo};
                    salFix = vertcat(fixItem.fix{idxOther});
                    for hf = 1:size(humanFix,1)
                        px = humanFix(hf,1); py = humanFix(hf,2);
                        if px > 0 && py > 0 && ...
                                px < size(fixItem.map,2) && py < size(fixItem.map,1)
                            humanMap(py,px) = humanMap(py,px) + 1;
                        end
                    end
                    for sf = 1:size(salFix,1)
                        px = salFix(sf,1); py = salFix(sf,2);
                        if px > 0 && py > 0 && ...
                                px < size(fixItem.map,2) && py < size(fixItem.map,1)
                            salMap(py,px) = salMap(py,px) + 1;
                        end
                    end
                    salMap = conv2(salMap,gauss,'same');
                    
                    %%% AUC part %%%%%%%%%%
                    % Eye fix vs shuffled eye fix of other img
                    [X Y] = find(humanMap > 0);
                    if isempty(X); continue; end
                    idx = 1:size(itsAllStiFiles,1); idx = idx(idx ~= i);
                    randIdx = randperm(length(idx));
                    idx = idx(randIdx(1:10)); % only take random 10 pics instead of all 1003
                    otherHumanMap = zeros(size(humanMap));
                    for zidx=idx
                        fixItemZ = get_NUSEF_fix_data(itsAllStiFiles(zidx).name);
                        otherHumanMap = otherHumanMap + imresize(im2double(fixItemZ.map),size(otherHumanMap),'nearest');
                    end
                    [XRest YRest] = find(otherHumanMap ~= 0);
                    
                    % Randomized pick
                    localHum = nan(length(X),1);
                    localRan = nan(length(X),100);
                    for k=1:length(X)
                        localHum(k,1) = salMap(X(k),Y(k)); % be aware
                        for kk=1:100
                            r = randi([1 length(XRest)],1);
                            localRan(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                        end
                        row = row + 1;
                    end
                    % scoring
                    [acz,R] = getauc(localHum, localRan, '');
                    tmpAC(subjNo) = mean(acz);
                end
                ac = tmpAC;
            else
                humanMap=im2double(fixItem.map);
                salMap = double(imresize(salMap,size(humanMap),'nearest'));
                %%% Gaussing
                mapSize = size(salMap);
                kSize = mapSize(2)*allGauss(gaussSize);
                salMap = imfilter(salMap, fspecial('gaussian', round([kSize, kSize]*4), kSize));
%                 salMap = conv2(salMap,gauss,'same');
                salMap = salMap - min(min(salMap));
                salMap = salMap / max(max(salMap));
                
                
                %%% AUC part %%%%%%%%%%
                [X Y] = find(humanMap > 0);
                idx = 1:size(itsAllStiFiles,1); idx = idx(idx ~= i);
                randIdx = randperm(length(idx));
                idx = idx(randIdx(1:10)); % only take random 10 pics instead of all
                otherHumanMap = zeros(size(humanMap));
                for zidx=idx
                    fixItem = get_NUSEF_fix_data(itsAllStiFiles(zidx).name);
                    otherHumanMap = otherHumanMap + imresize(im2double(fixItem.map),size(otherHumanMap),'nearest');
                end
                [XRest YRest] = find(otherHumanMap ~= 0);
                %             X = sacData(find(sacData(:,1)==imgIndex),2);
                %             Y = sacData(find(sacData(:,1)==imgIndex),3);
                %
                %             XRest = sacData(find(sacData(:,1)~=imgIndex),2);
                %             YRest = sacData(find(sacData(:,1)~=imgIndex),3);
                %
                localHum = nan(length(X),1);
                localRan = nan(length(X),100);
                for k=1:length(X)
                    %hum(row,1) = salMap(X(k),Y(k)); % be aware
                    localHum(k,1) = salMap(X(k),Y(k)); % be aware
                    for kk=1:100
                        r = randi([1 length(XRest)],1);
                        %ran(row,kk) = salMap(XRest(r),YRest(r)); % be aware
                        localRan(k,kk) = salMap(XRest(r),YRest(r)); % be aware
                    end
                    row = row + 1;
                end
                [ac,R] = getauc(localHum, localRan, '');
                
                
                %         [X Y] = find(humanMap ~= 0);
                %         NSSVector = zeros(1,size(X,1));
                %         for p=1:size(X,1)
                %             NSSVector(p) = salMap(X(p),Y(p));
                %         end
                %         r = mean(NSSVector);
                %r = corr2(salMap, humanMap);
                %             [ttr ttp] = corrcoef(salMap,humanMap);
                %         coefficients = [coefficients r];
                
            end
            
            save(['/lab/tmpir/u/dicky/nusef/cache/' myModelName '_' imageName '_' num2str(allGauss(gaussSize),'%1.2f') '.mat'],'-v7.3','ac')
        end
        coefficients(i) = mean(ac);
        %fprintf('.')
        %         if mod(i,100)==0; fprintf('%6ds | %s %.4d/%.4d\n',round(toc(myTic)),myModelName,i,length(list)); end
        fprintf('%6ds | %s %.4d/%.4d\n',round(toc(myTic)),myModelName,i,length(itsAllStiFiles));
        
        %         figure;
        %         imshow(salMap,[])
        %         pause
        clear salMap;
        
    end
    
    % correlations(count).mean = mean(coefficients(~isnan(coefficients)));
    % correlations(count).std = std(coefficients(~isnan(coefficients)));
    % correlations(count).coef = coefficients;
    % correlations(count).name = myModelName;
    
    addpath('/lab/raid/models/OutputResults')
    %     [ac,R] = getauc(hum, ran, '');
    %     correlations(count).mean = mean(coefficients);
    %     correlations(count).std = std(coefficients);
    %     correlations(count).coef = coefficients;
    %     correlations(count).name = myModelName;
    AUCValues(count).mean = mean(ac);
    AUCValues(count).std = std(ac);
    %     AUCValues(count).ttp = ttp;
    %     AUCValues(count).R = R;
    AUCValues(count).ac = coefficients;
    %     AUCValues(count).hum=hum;
    %     AUCValues(count).ran=ran;
    AUCValues(count).name = myModelName;
    AUCValues(count).coef = coefficients;
    %     count = count + 1;
    clear hum ran ac R;
    
    count = count + 1;
    
    fprintf('\n')
    
end
% error('balh')


myNum = ones(1,length(AUCValues)) * (length(itsAllStiFiles));


for ii=1:length(AUCValues)
    tmp = sum(~isnan(AUCValues(ii).coef));
    if tmp ~= myNum(ii); fprintf('Model %s has only %d out of %d\n', AUCValues(ii).name, tmp, myNum(ii)); end
    myNum(ii) = tmp;
    AUCValues(ii).mean = mean(AUCValues(ii).coef(~isnan(AUCValues(ii).coef)));
    AUCValues(ii).std = std(AUCValues(ii).coef(~isnan(AUCValues(ii).coef)));
end


xx = vertcat(AUCValues.mean);
yy = vertcat(AUCValues.std) ./ sqrt(myNum');
nn = {AUCValues.name};

    
nusefScores{1,gaussSize} = xx;
nusefScores{2,gaussSize} = yy;
nusefScores{3,gaussSize} = nn;
    
end

cd('/lab/raid/models/OutputResults');
save('./NewBottom-UpModel/Nov13/nusefLG.mat','-v7.3','nusefScores')
% save('aucAllScores.mat','-v7.3','bruceScores','juddScores','kootstraScores','nusefScores')


end
