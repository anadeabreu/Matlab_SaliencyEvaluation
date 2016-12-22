close all 

imgG = imread('C:\Users\Ana\Code\webvr-samples-master\panorama-viewer-master\Saliency\Results_20s_DBSCAN\GroupSaliencyMap_outdoor17.jpg'); % Ground Truth
FixG = imread('C:\Users\Ana\Code\webvr-samples-master\panorama-viewer-master\Saliency\Results_20s_DBSCAN\Fixations\GroupFixations_outdoor17.jpg'); % Ground Truth
%names = ['mergedSaliency_indoor3.jpg', 'indoor3_merged3.jpg']

%imgP = imread('C:\Users\Ana\Code\Images\QoMex\Salicon\mergedSaliency_outdoor17.jpg');
imgP = imread('C:\Users\Ana\Code\Images\QoMex\Salicon\outdoor17_merged3.jpg');

%imgP=rgb2gray(imgP);
imgG=rgb2gray(imgG);

[scoreNSS, meanAUC] = calcAUCscore (imgP, imgG);

FixG=rgb2gray(FixG);
BinFixG = imbinarize(FixG,0.1); % weird if T=0

[scoreNSS, meanNSS] = calcNSSscore(imgP,BinFixG);

meanAUC
meanNSS