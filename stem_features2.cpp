// NIST-developed software is provided by NIST as a public service. You may use,
// copy and distribute copies of the software in any medium, provided that you keep
// intact this entire notice. You may improve, modify and create derivative works of
//the software or any portion of the software, and you may copy and distribute such
// modifications or works. Modified works should carry a notice stating that you
// changed the software and should note the date and nature of any such change.
// Please explicitly acknowledge the National Institute of Standards and Technology
// as the source of the software.

// NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF
// ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT
// LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
// NON-INFRINGEMENT AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE
// OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS
// WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE
// USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS,
// ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

// You are solely responsible for determining the appropriateness of using and distributing
// the software and you assume all risks associated with its use, including but not limited
// to the risks and costs of program errors, compliance with applicable laws, damage to or
// loss of data, programs or equipment, and the unavailability or interruption of operation.
// This software is not intended to be used in any situation where a failure could cause risk
// of injury or damage to property. The software developed by NIST employees is not subject
// to copyright protection within the United States.

#include<iostream>
#include<fstream>
#include<math.h>
#include <tiffio.h>
#include <cstring>

using namespace std;

int findMean(uint16* pixeldata, uint16* mask, int n, int m) {
  double mean = 0.0;
  int cnt = 0;
  for (int i = 0; i < m*n; ++i) {
    if (mask[i] > 0) {
      mean += (double)(pixeldata[i]);
      ++cnt;
    }
  }
  if (cnt > 0){
    mean = mean/(double)(cnt);
  }

  double stdv = 0.0;
  for (int i = 0; i < m*n; ++i) {
    if (mask[i] > 0) {
      stdv += (double)(pixeldata[i]-mean)*(double)(pixeldata[i]-mean);
    }
  }
  if (cnt > 0){
    stdv = sqrt(stdv/(double)(cnt));
  }
  cout << "mean: " << mean << " standard deviation: " << stdv << endl;

  // mean of the whole image
  mean = 0.0;
  for (int i = 0; i < m*n; ++i) {
    mean += (double)(pixeldata[i]);
  }
  mean = mean/(double)(n*m);
  
  return mean;
}

double localstddev(uint16* pixel, uint16* mask, int width, int height) {
  int breakup = 10;
  float lmean, lsd;
  float smin = 10000.0;
  float smax = 0;
  double localstddev = 0.0;
  // need extra space to store pixel stddev values
  double *sddev = new double[width*height];
  int cc = 0;
  int cnt = 0;

  // collect sddev values for all pixels in the image
  for (int j = 0; j < height; ++j) {
    int firstj = j-breakup;
    int nextj = j+breakup;
    if (firstj < 0) firstj = 0;
    if (nextj > (height-1)) nextj = height-1;
    for (int k = 0; k < width; ++k) {
      cnt += 1;
      int firstk = k-breakup;
      int nextk = k+breakup;
      if (firstk < 0) firstk = 0;
      if (nextk > (width-1)) nextk = width-1;

      // collect these pixels
      cc = 0;
      lmean = 0.0;
      lsd = 0.0;
      int place = 0;
      for (int jj = firstj; jj < nextj; ++jj) {
	for (int ii = firstk; ii < nextk; ++ii) {
	  place = jj*width + ii;
	  lmean += (double)(pixel[place]);
	  ++cc;
	}
      }
      lmean = lmean/(double)(cc);
      sddev[j*width+k] = 0.0;
      for (int jj = firstj; jj < nextj; ++jj) {
	for (int ii = firstk; ii < nextk; ++ii) {
	  place = jj*width + ii;
	  sddev[j*width+k] += (lmean-(double)(pixel[place]))*(lmean-(double)(pixel[place]));
	}
      }
      sddev[j*width+k] = sqrt(sddev[j*width+k]/(double)(cc));
      if (sddev[j*width+k] < smin) smin = sddev[j*width+k];
      if (sddev[j*width+k] > smax) smax = sddev[j*width+k];
    }
  }

  lmean = 0.0;
  lsd = 0.0;

  // find the mean and stddev of the sddev background image (mask = 0)
  cc = 0;
  int cc2 = 0;
  for (int j = 0; j < height; ++j) {
    for (int k = 0; k < width; ++k) {
      if (mask[j*width+k] < 1) {
	lmean += sddev[j*width+k];
	cc += 1;
      }
      else {
	cc2 += 1;
      }
    }
  }
  lmean = lmean/(double)(cc);
  for (int j = 0; j < height; ++j) {
    for (int k = 0; k < width; ++k) {
      if (mask[j*width+k] < 1) {
	lsd += (lmean-sddev[j*width+k])*(lmean-sddev[j*width+k]);
      }
    }
  }
  lsd = sqrt(lsd/(float)(cc));

  // find 3 sd above the background mean
  double bkmean3 = lmean + 3.0*lsd;
  //cout << "3 sds above the background mean sddev: " << bkmean3 << endl;

  // find fraction of cell pixels above bkmean3;
  cc = 0;
  cc2 = 0;
  for (int j = 0; j < height; ++j) {
    for (int k = 0; k < width; ++k) {
      if (mask[j*width+k] > 0) {
	if (sddev[j*width+k] > bkmean3) {
	  cc += 1;
	}
	else {
	  cc2 += 1;
	}
      }
    }
  }

  localstddev = 1.0;
  if (cc2 > 0)
    localstddev = double(cc)/double(cc+cc2);
  //cout << "final: " << cc << " " << cc2 << " " << localstddev << endl;
  delete [] sddev;
  return localstddev;
}


double fgbg(uint16* pixel, uint16* mask, int width, int height) {
  double fgbgval = 0.0;
  // need extra space to for perimeter dilation
  int *grid = new int[width*height];
  int *grid2 = new int[width*height];
  for (int j = 0; j < width*height; ++j) {
    grid[j] = 0;
    grid2[j] = 0;
  }

  // put the perimeter into grid1
  for (int x = 0; x < height; x++) {
    for (int y = 0; y < width; y++) {
      if (mask[x*width+y] > 0) {
	if (x == 0 || y == 0 || y == width-1 || x == height-1) {
	  //some neighbors will be out
	  grid[x*width+y] = 1;
	}
	// look for an edge (mask = 0)
	else if ((mask[x*width+y-1] < 1) || (mask[x*width+y+1] < 1) || (mask[x*width+y-width] < 1) || (mask[x*width+y+width] < 1)) {
	  grid[x*width+y] = 1;
	}
	else if ((mask[x*width+y-1-width] < 1) || (mask[x*width+y+1-width] < 1) || (mask[x*width+y-1+width] < 1) || (mask[x*width+y+1+width] < 1)) {
	  grid[x*width+y] = 1;
	}
      }
    }
  }

  // find the area 50 pixels inside (grid=2) and outside (grid=3) the perimeter
  for (int times = 0; times < 50; ++times) {
    // dilate: first copy to grid2
    for (int j = 0; j < width*height; ++j) {
      grid2[j] = grid[j];
    }

    // dilate
    for (int x = 0; x < height; x++) {
      for (int y = 0; y < width; y++) {
	if (grid[x*width+y] > 0) {
	  grid2[x*width+y-1] = 2;
	  grid2[x*width+y+1] = 2;
	  grid2[x*width+y-width] = 2;
	  grid2[x*width+y+width] = 2;
	}
      }
    }

    // copy back to grid
    for (int j = 0; j < width*height; ++j) {
      grid[j] = grid2[j];
    }    
  }

  // mark inside or outside
  for (int x = 0; x < height; x++) {
    for (int y = 0; y < width; y++) {
      if (grid[x*width+y] > 1) {
	if (mask[x*width+y] > 0)
	  grid[x*width+y] = 2;
	else
	  grid[x*width+y] = 3;
      }
    }
  }
  
  // find the ratio of the entropies of areas 1 and 2
  // find the max pixel and get histograms
  int pmax = 0;
  for (int j = 0; j < width*height; ++j) {
    if (pixel[j] > pmax) pmax = pixel[j];
  }

  double *hist1 = new double[pmax];
  double *hist2 = new double[pmax];
  for (int j = 0; j < pmax; ++j) {
    hist1[j] = 0.0;
    hist2[j] = 0.0;
  }
  
  uint16 val;
  int cnt2 = 0;
  int cnt3 = 0;
  for (int x = 0; x < height; x++) {
    for (int y = 0; y < width; y++) {
      if (grid[x*width+y] == 2) {
	val = pixel[x*width+y];
	hist1[val] += 1.0;
	cnt2 += 1;
      }
      if (grid[x*width+y] == 3) {
	val = pixel[x*width+y];
	hist2[val] += 1.0;
	cnt3 += 1;
      }
    }
  }

  // normalize histograms  
  for (int p = 0; p < pmax; ++p) {
    hist1[p] = hist1[p]/(double)(cnt2);
    hist2[p] = hist2[p]/(double)(cnt3);
  }

  // get entropy of the inside 50 pixels
  double entropy1 = 0.0;
  double entropy2 = 0.0;
  double maxlog = -log(1.0/(double)(pmax));
  double rangelog = maxlog/100.0;
  
  double eee;
  for (int i = 0; i < pmax; ++i) {
    if (hist1[i] > 0) {
      eee = hist1[i] * log(hist1[i]);
      entropy1 += 100.0 - (maxlog - eee)/rangelog;
    }
    if (hist2[i] > 0) {
      eee = hist2[i] * log(hist2[i]);
      entropy2 += 100.0 - (maxlog - eee)/rangelog;
    }
  }

  fgbgval = entropy1;
  if (fabs(entropy2) > 0.0000001) fgbgval = entropy1/entropy2;

  delete [] grid;
  delete [] grid2;
  delete [] hist1;
  delete [] hist2;
  return fgbgval;
}

int main(int argc, char *argv[]) {
/// \param argv arguments list
///             argv[0] program name
///             argv[1] path mask
///             argv[2] path image

  // read in the mask
  TIFF *tif=TIFFOpen(argv[1], "r");
  uint16 *mask;
  int width;
  int height;
  uint16 row, col;
  uint16 minval = 10000;
  uint16 maxval = 0;
  if (tif) {
      TIFFGetField(tif,TIFFTAG_IMAGEWIDTH, &width);
      TIFFGetField(tif,TIFFTAG_IMAGELENGTH, &height);

      mask = new uint16[width*height];
      for (row = 0; row < height; row++) {
	TIFFReadScanline(tif, &mask[row*width], row);
	for (col = 0; col < width; col++) {
	  if (mask[row*width+col] < minval) minval = mask[row*width+col];
	  if (mask[row*width+col] > maxval) maxval = mask[row*width+col];
	}
      }
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
  }

  // avg intensities
  findMean(pix,mask,width,height);

  // find localstddev

  // get a new mask with just 0 and 1's
  for (int j = 0; j < width*height; ++j) {
    if (mask[j] > 0)
      mask[j] = 1;
  }
  double locsd = localstddev(pix,mask,width,height);
  cout << "found localstddev: " << locsd << endl;

  // find fgbg
  double fgbgval = fgbg(pix,mask,width,height);
  cout << "found fbgb: " << fgbgval << endl;

  delete(pix);
  delete(mask);
  delete(buff16);

  return 0;
}
