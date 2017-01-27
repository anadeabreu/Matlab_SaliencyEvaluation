close all 

%image_names = {'BellagioHotelLobby', 'CathedraledeBordeaux'}

image_names = {'BellagioHotelLobby', 'CathedraledeBordeaux', 'theBeach', 'OverParis', 'indoor3', 'indoor4', 'indoor5', 'indoor6', 'indoor7', 'outdoor1', 'outdoor2', 'outdoor4', 'outdoor5', 'outdoor6', 'outdoor8', 'outdoor9', 'outdoor15', 'outdoor16', 'outdoor17', 'outdoor18',  'The_Porch'}
fprB_images = [];
tprB_images = [];
scores_imagesB = [];
scores_imagesJ = [] ;
scores_imagesS = [] ;
scoresNSS_images = [] ;

for i=1: size(image_names,2)
    


    %imgG2 = imread('C:\Users\Ana\Code\webvr-samples-master\panorama-viewer-master\Saliency\10s\bwGroupSaliency_The_Porch.jpg'); % Ground Truth
    FixG = imread(['C:\Users\Ana\Code\webvr-samples-master\panorama-viewer-master\Saliency\10s\GroupFix_' image_names{1,i} '_0.jpg']); % Ground Truth
    otherMap = imread('C:\Users\Ana\Code\webvr-samples-master\panorama-viewer-master\Saliency\10s\GroupFix_indoor6_0.jpg');

    %imgP = imread('C:\Users\Ana\Code\Images\QoMex\Salicon\mergedSaliency_outdoor5.jpg');
    imgP = imread(['C:\Users\Ana\Code\Images\QoMex\Salicon\mergedSaliency_' image_names{1,i} '.jpg']);

    %imgP = rgb2gray(imgP);
    imgP = mat2gray(imgP);
 
 
    FixG = rgb2gray(FixG);
    FixG = mat2gray(FixG);
 
     BinFixG = imbinarize(FixG,0.1);
 
    otherMap = rgb2gray(otherMap);
    otherMap = mat2gray(otherMap);
 
    BinOtherMap = imbinarize(otherMap,0.1);

% imgG1=rgb2gray(imgG1);
% imgG2=rgb2gray(imgG2);
% BinImgG1 = imbinarize(imgG1,0.1);
% BinImgG2 = imbinarize(imgG2,0.1);
% 
%  size(imgP)
%  size(FixG)

    [scoreJ,tpJ,fpJ,allthreshes] = AUC_Judd (imgP, BinFixG, 0, 1);


    [scoreB,tpB,fpB] = AUC_Borji (imgP, BinFixG, 100, 0.01, 1);
% length(tp)
% allOneString = sprintf('%.4f,' , tp');
% allOneString = allOneString(1:end-1)
% scoreB

    [scoreS,tp,fp] = AUC_shuffled(imgP, BinFixG, BinOtherMap);
% scoreS

    [scoreNSS] = NSS(imgP,BinFixG)
    
    length(fpB)
    
    while (length(fpB)<102)
        fpB = [fpB; 1]
        tpB = [tpB; 1]
    end
    
    % Using Borji for the plot
    fpB'
    tpB'
    fprB_images = [fprB_images; fpB(1:102)'];
    tprB_images = [tprB_images; tpB(1:102)'];
    scores_imagesB = [scores_imagesB scoreB] ;
    scores_imagesJ = [scores_imagesJ scoreJ] ;
    scores_imagesS = [scores_imagesS scoreS] ;
    scoresNSS_images = [scoresNSS_images scoreNSS] ;
end

 fprB_images = mean(fprB_images)
 tprB_images = mean(tprB_images)
% scores_imagesB
% scores_imagesJ
% scores_imagesS
% scoresNSS_images
 
% meanAUC
% meanNSS