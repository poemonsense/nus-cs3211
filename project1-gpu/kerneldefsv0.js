
var gpu = new GPU();

function sqr(x) {
    return x*x;
}

function dist(x1,y1,x2,y2) {
    return Math.sqrt( sqr(x2-x1)+sqr(y2-y1) );
}

gpu.addFunction(sqr);
gpu.addFunction(dist);

// don't what this is about
// assume it does nothing
function makeImg(mode) {
    var opt = {
        dimensions: [800, 600, 4],
        graphical: false,
        outputToTexture: true,
        mode: mode
    };
    var y = gpu.createKernel(function(img) {
        return img[this.thread.z][this.thread.y][this.thread.x];
    }, opt);
    return y;
}

var toimg = gpu.createKernel(function(A) {
    this.color(A[0][this.thread.y][this.thread.x],A[1][this.thread.y][this.thread.x],A[2][this.thread.y][this.thread.x]);
}).dimensions([800, 600]).graphical(true);

// demo filter provided by Handsome Hugh
function makeFilter(mode) {
    if (mode === "cpu")
        var dim = [800, 600];
    else
        var dim = [800, 600, 3];
    var opt = {
        dimensions: dim,
        graphical: false,
        outputToTexture: true,
        mode: mode
    };
    var filt = gpu.createKernel(function(A) {
        if (this.thread.y > 0 && this.thread.y < 599 && this.thread.x < 799 && this.thread.x > 0) {
            // group the same terms and produce a simplified form
            var c = A[this.thread.y][this.thread.x+1] + 
                    A[this.thread.y+1][this.thread.x+1] +
                    A[this.thread.y+1][this.thread.x] - 
                    A[this.thread.y-1][this.thread.x-1] -
                    A[this.thread.y][this.thread.x-1] -
                    A[this.thread.y-1][this.thread.x];
            return 2*c+0.5;
        } else {
            return A[this.thread.y][this.thread.x];
        }
    },opt);
    // it seems that in cpu mode, this lambda function works pretty good
    if (mode == "cpu")
        return x => [filt(x[0]), filt(x[1]), filt(x[2])];
    else
        return filt;
}

// change an RGB image into grayscale (returns just one channel)
function gray_filter(mode, dim) {
    var opt = {
        dimensions: dim,
        graphical: false,
        outputToTexture: true,
        mode: mode
    };
    // the Y'UV and Y'IQ models used by PAL and NTSC
    // ref: https://en.wikipedia.org/wiki/Grayscale
    var mask = gpu.createKernel(function(A) {
        return 0.299*A[0][this.thread.y][this.thread.x]+
               0.587*A[1][this.thread.y][this.thread.x]+
               0.114*A[2][this.thread.y][this.thread.x];
    },opt);
    return mask;
}

// change an RGB image into grayscale (return 800*600*3 RGB image as well)
function gray_filter_3d(mode) {
    if (mode === "cpu")
        return x => {
            var y = gray_filter(mode, [800,600])(x);
            return [y,y,y];
        };
    else
        return gray_filter(mode, [800,600,3]);
}

// gaussian blur for any channel with image size of 800*600
// parameters for gaussian blurring: size = 3*3, sigma = 2
// use the sample matrix instead of the brute-force computation
function gaussian_blur_1d(mode, dim) {
    return gpu.createKernel(function(A) {
        if (this.thread.y > 0 && this.thread.y < this.constants.height-1 &&
                this.thread.x < this.constants.width-1 && this.thread.x > 0) {
            return  A[this.thread.y-1][this.thread.x-1]*0.10187 +
                    A[this.thread.y-1][this.thread.x]*0.11543 +
                    A[this.thread.y-1][this.thread.x+1]*0.10187 +
                    A[this.thread.y][this.thread.x-1]*0.11543 +
                    A[this.thread.y][this.thread.x]*0.13080 +
                    A[this.thread.y][this.thread.x+1]*0.11543 +
                    A[this.thread.y+1][this.thread.x-1]*0.10187 +
                    A[this.thread.y+1][this.thread.x]*0.11543 +
                    A[this.thread.y+1][this.thread.x+1]*0.10187;
        } else {
            return A[this.thread.y][this.thread.x];
        }
    }).dimensions(dim).outputToTexture(true).mode(mode).constants({width:dim[0], height:dim[1]});
}

// perform gaussian blur on the 800*600*3 RGB image
function gaussian_blur(mode) {
    if (mode === "cpu") {
        var filt = gaussian_blur_1d(mode, [800,600]);
        return x => [filt(x[0]), filt(x[1]), filt(x[2])];
    }
    else {
        return gpu.createKernel(function(A) {
            if (this.thread.y > 0 && this.thread.y < 599 &&
                    this.thread.x < 799 && this.thread.x > 0) {
                return  A[this.thread.z][this.thread.y-1][this.thread.x-1]*0.10187 +
                        A[this.thread.z][this.thread.y-1][this.thread.x]*0.11543 +
                        A[this.thread.z][this.thread.y-1][this.thread.x+1]*0.10187 +
                        A[this.thread.z][this.thread.y][this.thread.x-1]*0.11543 +
                        A[this.thread.z][this.thread.y][this.thread.x]*0.13080 +
                        A[this.thread.z][this.thread.y][this.thread.x+1]*0.11543 +
                        A[this.thread.z][this.thread.y+1][this.thread.x-1]*0.10187 +
                        A[this.thread.z][this.thread.y+1][this.thread.x]*0.11543 +
                        A[this.thread.z][this.thread.y+1][this.thread.x+1]*0.10187;
            } else {
                return A[this.thread.z][this.thread.y][this.thread.x];
            }
        }).dimensions([800,600,3]).outputToTexture(true).mode(mode);
    }
}

// change the 800*600 image into a 200*150 image
function reduce_size(mode) {
    return gpu.createKernel(function(A) {
        return A[this.thread.y*4][this.thread.x*4];
    }).dimensions([200, 150]).outputToTexture(true).mode(mode);
}

// enlarge the 200*150 grayscale image to 800*600*3 RGB image
function enlarge(mode) {
    return gpu.createKernel(function(A) {
        return A[Math.floor(this.thread.y/4)][Math.floor(this.thread.x/4)];
    }).dimensions([800, 600, 3]).outputToTexture(true).mode(mode);
}

var sin_val = [];
var cos_val = [];
for (var i = 0; i < 36; i++) {
    var rad = (5*i-90)*Math.PI/180;
    sin_val[i] = Math.sin(rad);
    cos_val[i] = Math.cos(rad);
}

// soble operator for edge detection
function sobel(mode,dim) {
    return gpu.createKernel(function(A) {
        // ref: https://en.wikipedia.org/wiki/Sobel_operator
        if (this.thread.y > 0 && this.thread.y < this.constants.height-1 &&
                this.thread.x < this.constants.width && this.thread.x > 0) {
            var gx = A[this.thread.y+1][this.thread.x-1] + A[this.thread.y+1][this.thread.x+1] +
                     2*(A[this.thread.y+1][this.thread.x] - A[this.thread.y-1][this.thread.x]) -
                     (A[this.thread.y-1][this.thread.x-1] + A[this.thread.y-1][this.thread.x+1]);
            var gy = A[this.thread.y-1][this.thread.x+1]-A[this.thread.y-1][this.thread.x-1]+
                     2*(A[this.thread.y][this.thread.x+1]-A[this.thread.y][this.thread.x-1])+
                     A[this.thread.y+1][this.thread.x+1]-A[this.thread.y+1][this.thread.x-1];
            var g = Math.sqrt(gx*gx + gy*gy);
            return g;
        }
        return 0;
    }).dimensions(dim).outputToTexture(true).mode(mode).constants({width:dim[0],height:dim[1]});
}

// first mask for line detector (use hough transform)
// loop over all white pixels and implement the accumulator
// it seems to take a lot of time to compile via gpu.rocks
// probably due to the various loop bound 
function line_detector_mask1(mode) {
    return gpu.createKernel(function(A, cos_val, sin_val) {
        var count = 0;
        if (this.thread.x != 0) {
            for (var i = 0; i < this.constants.height; i++) {
                var resi = this.thread.y - 250 - i * sin_val[this.thread.x];
                var j = Math.floor(resi / cos_val[this.thread.x]);
                if (A[i][j] > 0.6)
                    count += 1;
            }
        }
        // horizontal line
        else if (this.thread.y >= 250 && this.thread.y < 400) {
            var i = this.thread.y - 250;
            for (var j = 0; j < this.constants.width; j++) {
                if (A[i][j] > 0.6)
                    count += 1;
            }
        }
        return count;
    }).dimensions([36, 500]).outputToTexture(true).mode(mode).constants({width:200,height:150});
}

// second mask for line detector (use hough transform)
// find all the (r, theta) with satisfiable accumulator result
// and make all the pixels on that line white
function line_detector_mask2(mode) {
    return gpu.createKernel(function(A, map, cos_val, sin_val) {
        for (var j = 0; j < this.constants.size; j++) {
            // compute r instead of looping over all possibilties
            var r = Math.floor(this.thread.x * cos_val[j] + this.thread.y * sin_val[j]) + 250;
            if (map[r][j] > 60) {
                return 1;
            }
        }
        return A[this.thread.y][this.thread.x];
    }).dimensions([200,150]).outputToTexture(true).mode(mode).constants({size:36});
}

// pixelate the image
// i.e. replace all 4*4 square with the top left pixel
function pixelated(mode) {
    return gpu.createKernel(function(A) {
        var y = this.thread.y - (this.thread.y % 4);
        var x = this.thread.x - (this.thread.x % 4);
        return A[this.thread.z][y][x];
    }).dimensions([800, 600, 3]).outputToTexture(true).mode(mode);
}

// sharpening
// http://northstar-www.dartmouth.edu/doc/idl/html_6.2/Sharpening_an_Image.html
function sharpen(mode) {
    return gpu.createKernel(function(A) {
        if (this.thread.y > 0 && this.thread.y < 599 &&
                this.thread.x < 799 && this.thread.x > 0) {
            var c = 9*A[this.thread.z][this.thread.y][this.thread.x]-
                    (A[this.thread.z][this.thread.y-1][this.thread.x-1] +
                    A[this.thread.z][this.thread.y-1][this.thread.x] +
                    A[this.thread.z][this.thread.y-1][this.thread.x+1] +
                    A[this.thread.z][this.thread.y][this.thread.x-1] +
                    A[this.thread.z][this.thread.y][this.thread.x+1] +
                    A[this.thread.z][this.thread.y+1][this.thread.x-1] +
                    A[this.thread.z][this.thread.y+1][this.thread.x] +
                    A[this.thread.z][this.thread.y+1][this.thread.x+1]);
            return A[this.thread.z][this.thread.y][this.thread.x] + c / 9;
        } else {
            return A[this.thread.z][this.thread.y][this.thread.x];
        }
    }).dimensions([800,600,3]).outputToTexture(true).mode(mode);
}