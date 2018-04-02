mkdir bmpresults
cd ppmresults
find . -name '*.ppm' |  while read line; do basename $line .ppm; done | while read line; do echo "Processing $line.ppm"; ppmtobmp "$line.ppm" > "../bmpresults/$line.bmp"; done
