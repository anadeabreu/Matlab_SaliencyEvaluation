close all 

imgG1 = imread('C:\Users\Ana\Code\webvr-samples-master\panorama-viewer-master\Saliency\10s\bwGroupSaliency_BellagioHotelLobby.jpg'); % Ground Truth
imgG2 = imread('C:\Users\Ana\Code\webvr-samples-master\panorama-viewer-master\Saliency\10s\bwGroupSaliency_The_Porch.jpg'); % Ground Truth
%FixG = imread('C:\Users\Ana\Code\webvr-samples-master\panorama-viewer-master\Saliency\10s\Fixations\GroupFixations_outdoor17.jpg'); % Ground Truth
%names = ['mergedSaliency_indoor3.jpg', 'indoor3_merged3.jpg']

%imgP = imread('C:\Users\Ana\Code\Images\QoMex\Salicon\mergedSaliency_outdoor17.jpg');
imgP = imread('C:\Users\Ana\Code\Images\QoMex\Salicon\mergedSaliency_The_Porch.jpg');

imgP=rgb2gray(imgP);
imgG1=rgb2gray(imgG1);
imgG2=rgb2gray(imgG2);

BinImgG1 = imbinarize(imgG1,0.1);
BinImgG2 = imbinarize(imgG2,0.1);
imshow(BinImgG1)

[score, meanAUC] = calcAUCscore (imgP, BinImgG2, BinImgG1);

% FixG=rgb2gray(FixG);
% BinFixG = imbinarize(FixG,0.1); % weird if T=0
% 
% [scoreNSS, meanNSS] = calcNSSscore(imgP,BinFixG);
% 
% meanAUC
% meanNSS