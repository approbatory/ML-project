# Noise correlations in neural ensemble activity limit the accuracy of hippocampal spatial representations

The code used for the paper by Omer Hazon, Victor Minces, David Tomàs, Surya Ganguli, Mark Schnitzer, and Pablo Jercog.

## Selected functions implementing our methods

- Code used for detecting events in calcium traces: `utils/hyperdetect.m`
- Code used for computing the position signal and noise along the signal direction: `DecodeTensor.m` function `adjacent_metrics`
- Code used for computing the signal, noise, and SNR along the noise covariance eigenvectors: `@Cloud/eigen_snr_crossval.m`
- Code used for computing the loadings (cos angle) between the noise covariance eigenvectors and the signal direction, plus related quantities:`@Cloud/Cloud.m` function `Cloud`
- Code used for computing the normalized signal variance: `@Org/load_definitions.m` lambda function on line 43
