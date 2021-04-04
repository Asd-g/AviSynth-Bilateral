# Description

Bilateral filter, can be used to perform spatial de-noise, spatial smoothing while preserving edges.

Larger spatial sigma results in larger smoothing radius, while smaller range sigma preserves edges better.

Now there're 2 different algorithms implemented in this function, algorithm=1 is suitable for large sigmaS and large sigmaR, and algorithm=2 is suitable for small sigmaS and small sigmaR. By default, algorithm=0 will choose the appropriate algorithm to perform the filtering.

If clip ref is specified, it will perform joint/cross Bilateral filtering, which means clip ref is used to determine the range weight of Bilateral filter.

By default, this function will only process Y plane for YUV format, and process all the planes for other formats. It is recommended to use Y as ref clip for chroma plane processing.

This is [a port of the VapourSynth plugin Bilateral](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-Bilateral).

# Usage

```
Bilateral (clip input, clip "ref", float "sigmaSY", float "sigmaSU", float "sigmaSV", float "sigmaRY", float "sigmaRU", float "sigmaRV", int "algorithmY", int "algorithmU", int "algorithmV", int "PBFICnumY", int "PBFICnumU", int "PBFICnumV", int "y", int "u", int "v")
```

## Parameters:

- input\
    A clip to process.\
    It must be in 8..16-bit planar format.
    
- ref\
    Reference clip to calculate range weight.\
    Specify it if you want to perform joint/cross Bilateral filter.\
    Must have the same properties as input.\
    Default: input.
    
- sigmaSY, sigmaSU, sigmaSV\
    Sigma of Gaussian function to calculate spatial weight.\
    The scale of this parameter is equivalent to pixel distance.\
    Larger values results in larger filtering radius as well as stronger smoothing.\
    Must be non-negative.\
    algorithmX = 1: It is of constant processing time regardless of sigmaS, while in small sigmaS the smoothing effect is stronger compared to Bilateral filter prototype.\
    algorithmX = 2: It will be slower as sigmaS increases, and for large sigmaS the approximation will be bad compared to Bilateral filter prototype.\
    Default: sigmaSY = 3.0; sigmaSU/sigmaSV depending on the video chroma subsampling.

- sigmaRY, sigmaRU, sigmaRV\
    Sigma of Gaussian function to calculate range weight.\
    The scale of this parameter is the same as pixel value ranging in `[0,1]`.\
    Smaller sigmaR preserves edges better, may also leads to weaker smoothing.\
    Must be non-negative.\
    algorithmX = 1: As sigmaR decreases, the approximation of this algorithm gets worse, so more PBFICs should be used to produce satisfying result. If PBFICnum is not assigned, the number of PBFICs used will be set according to sigmaR.\
    algorithmX = 2: It is of constant processing time regardless of sigmaR, while for large sigmaR the approximation will be bad compared to Bilateral filter prototype.\
    Default: sigmaRY = sigmaRU = sigmaRV = 0.02.
    
- algorithmY, algorithmU, algorithmV\
    0: Automatically determine the algorithm according to sigmaS, sigmaR and PBFICnum.\
    1: O(1) Bilateral filter uses quantized PBFICs.\
    2: Bilateral filter with truncated spatial window and sub-sampling. O(sigmaS^2).\
    Default: algorithmY = algorithmU = algorithmV = 1.
    
- PBFICnumY, PBFICnumU, PBFICnumV\
    Number of PBFICs used in algorithm=1.\
    Must be between 2..256.\
    Default: 4 when sigmaRX >= 0.08. It will increase as sigmaR decreases, up to 32. For chroma plane default value will be odd to better preserve neutral value of chromiance.
    
- y, u, v\
    Planes to process.\
    1: Return garbage.\
    2: Copy plane.\
    3: Process plane. Always process planes when clip is RGB.\
    Default: y = 3; u = v = 1.
    
# Building

## Windows

Use solution files.

## Linux

### Requirements

- Git
- C++17 compiler
- CMake >= 3.16

```
git clone https://github.com/Asd-g/AviSynth-Bilateral && \
cd AviSynth-Bilateral && \
mkdir build && \
cd build && \
cmake .. && \
make -j$(nproc) && \
sudo make install
```
