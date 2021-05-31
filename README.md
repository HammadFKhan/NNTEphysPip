# NNTEphysPip

![gitsha-01](https://user-images.githubusercontent.com/56140216/117356034-50877180-ae81-11eb-8e78-e014f26c2eab.png)

Electrophysiology Pipeline

### Image and data processing pipeline for eletrophysiological recordings

The main focus of this pipline is to create and combine multiple repositories to comprehensivley analyze eletrophyiological data either from shanks or microelectrode arrays. This pipeline currentley consists of clusterless analysis. More elaborate explainations are found later in this document. 

## Installation

### Requirements

- MATLAB 2016     
- Signal Processing Toolbox

Please download the full repository into your working MATLAB folder. You must also download the reposity of [buzcode](https://github.com/buzsakilab/buzcode) for ROI extraction and motion correction of the image stack. You will need to add the [buzcode](https://github.com/buzsakilab/buzcode) folder to the NNTCaPip MATLAB path.

If your computer is CUDA compatible, I highly recommend you also download [Kilosort3](https://github.com/MouseLand/Kilosort) to analyize spiketrains using clustered sorting. 


