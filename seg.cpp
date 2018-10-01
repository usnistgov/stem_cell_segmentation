#include<iostream>
#include<fstream>
#include<math.h>
#include <tiffio.h>
#include <cstring>
using namespace std;
const float pi = 3.14159;

#include "segment.cpp"

int main(int argc, char *argv[]) {

  // read in the mask
  TIFF *tif=TIFFOpen(argv[1], "r");
  //int8_t *mask;
  uint16 *mask;
  int width;
  int height;
  uint16 row, col;
  uint16 minval = 10000;
  uint16 maxval = 0;
  if (tif) {
      TIFFGetField(tif,TIFFTAG_IMAGEWIDTH, &width);
      TIFFGetField(tif,TIFFTAG_IMAGELENGTH, &height);
      cout << width << " " << height << endl;

      mask = new uint16[width*height];
      for (row = 0; row < height; row++) {
	TIFFReadScanline(tif, &mask[row*width], row);
	for (col = 0; col < width; col++) {
	  if (mask[row*width+col] < minval) minval = mask[row*width+col];
	  if (mask[row*width+col] > maxval) maxval = mask[row*width+col];
	}
      }
      cout << minval << " " << maxval << endl;
  }

  uint16 *pix;
  tif=TIFFOpen(argv[2], "r");
  minval = 10000;
  maxval = 0;
  char *buff16 = new char[sizeof(uint16)*width];
  if (tif) {
      TIFFGetField(tif,TIFFTAG_IMAGEWIDTH, &width);
      TIFFGetField(tif,TIFFTAG_IMAGELENGTH, &height);

      pix = new uint16[width*height];
      for (row = 0; row < height; row++) {
	TIFFReadScanline(tif, buff16, row);
	for (col = 0; col < width; col++) {
	  memcpy(&pix[row*width+col], &buff16[2*col], sizeof(uint16));
	  if (pix[row*width+col] < minval) minval = pix[row*width+col];
	  if (pix[row*width+col] > maxval) maxval = pix[row*width+col];
	}
      }
      cout << minval << " " << maxval << endl;
  }

  segment *seg = new segment(pix,width,height);
  seg->findmask();
  cout << "found mask" << endl;

  delete(pix);
  delete(mask);
  delete(buff16);

  return 0;
}
