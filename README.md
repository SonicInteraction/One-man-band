# One-man-band
## Audio analysis, feature extraction, pitch estimation and automated time stretching and pitch shifting in matlab

__**Below is a list of files and a precis of their functionality**__

*setup.m*
- load files
- determine lengths
- load all into array
- audition
- load any previous time warped files
- mix (normalise)
- plot a waterfall spectrogram

*pre-processing.m*
- band pass filter
- compression
- amplitude plot for comparison

*pitch-detection.m*
- autocorrelation
- some pre-processing
- ***Needs tidying up***

*amplitude-frequency-analysis-time-warp-and save results.m*
- Calculate FFT
- plot
    - Amplitude
    - Spectrogram
    - dissimilarity matrix
    - cumulative dissimilarity matrix
    - visualisation of time warping path
 - apply time warp
 - amplitude plot of result and play results
