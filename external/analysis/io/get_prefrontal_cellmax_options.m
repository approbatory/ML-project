function options = get_prefrontal_cellmax_options()

options.movieDatasetName = '/Data/Images';
options.useParallel = 0;
options.useImageCorrEventDetect = 0;

options.CELLMaxoptions.inputSizeManual = 0;
options.CELLMaxoptions.gridWidth = 10;
options.CELLMaxoptions.gridSpacing = 18;
options.CELLMaxoptions.maxSqSize = 120;
options.CELLMaxoptions.sqOverlap = 20;
options.CELLMaxoptions.selectRandomFrames = 1;
options.CELLMaxoptions.numFramesRandom = 4000;
options.CELLMaxoptions.minIters = 52;
options.CELLMaxoptions.initMethod = 'ica';
