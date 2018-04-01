
var selection = 0;
var filtering = 0;
var counter = 0;

function change( el ) {
    if ( el.value === "Using CPU" ) {
        selection = 1;
        el.value = "Using GPU";
    } else {
        selection = 0;
        el.value = "Using CPU";
    }
}

function changeFilter( el ) {
    if ( el.value === "Filtering" ) {
        filtering = 1;
        el.value = "Gauss";
    } else if (el.value === "Gauss") {
        filtering = 2;
        el.value = "Pixelate";
    } else if (el.value === "Pixelate") {
        filtering = 3;
        el.value = "Sharpen"
    } else if (el.value === "Sharpen") {
        filtering = 4;
        el.value = "Edge Detector"
    } else if (el.value === "Edge Detector") {
        filtering = 5;
        el.value = "Line Detector"
    } else if (el.value === "Line Detector") {
        filtering = 6;
        el.value = "No Filter"
    } else {
        filtering = 0;
        el.value = "Filtering";
    }
}



