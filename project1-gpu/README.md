# Image Processing: Implementation and Analysis

Develop a small set of original image processing functions, which run both on a GPU and CPU
platform, in a browser. 

Doing this (only) in a subset of Javascript, using `GPU.js` (the original GPU accelerated JavaScript http://gpu.rocks/).

## Implementation

- image resizing filters
- a grayscale filter based on the Y'UV and Y'IQ models $Y = 0.299R + 0.587G + 0.114B$
- a pixelation filter
- an image sharpening filter based on unsharp masking (USM) using the kernel
- gaussian blurring with a $3\times 3$ Gaussian kernel (with $\sigma =2$ )
- sobel operator, which creates an image emphasising edges and is used for edge detection
- line detector based on the Hough transform of $r=x\cos\theta+y\sin\theta$ 