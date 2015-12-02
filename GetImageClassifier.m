function [ Classifier ] = GetImageClassifier( images, labels )
%% MakeImageClassifier
%   Creates a multi-label image classifier given the training set train_X
%   and labels train_Y
%   Input:
%       train_X - an n-by-m-by-p matrix containing p greyscale images that 
%       are all n by m pixels.
%       train_Y - a vector of length p containing the labels of each image.

%% Parameters
% HogCellSize -> Lower values encodes more information (higher accuracy)
HogCellSize = [2, 2];
p = size(images, 3);

%% Setup
% Get the HOG feature length from the first image.
firstImage = images(:,:,1);
featureLen = length(extractHOGFeatures(firstImage, 'CellSize', HogCellSize));
% Initialize the features array (training dataset)
features = zeros(p, featureLen, 'single');

%% Extract the HOG features for each image
for i = 1 : p
    image = images(:,:,i);
    
    grayscale = mat2gray(image, [0, max(max(image))]);
    features(i, :) = extractHOGFeatures(grayscale, 'CellSize', HogCellSize);
end

%% Train the classifier
Classifier = fitcecoc(features, labels);

end

