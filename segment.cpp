#include "segment.h"
#include "extras.cpp"

void segment::findmask() {

  // get the mean and sdv of the image
  meanv = 0.0;
  int cnt = 0;
  for (int i = 0; i < width*height; ++i) {
    meanv += (double)(pixel[i]);
    ++cnt;
  }
  if (cnt > 0){
    meanv = meanv/(double)(cnt);
  }

  sdv = 0.0;
  for (int i = 0; i < width*height; ++i) {
    sdv += (double)(pixel[i]-meanv)*(double)(pixel[i]-meanv);
  }
  if (cnt > 0){
    sdv = sqrt(sdv/(double)(cnt));
  }
  cout << "mean: " << meanv << " standard deviation: " << sdv << endl;

  // change zeros so they do not mess up segmentation
  for (int i = 0; i < width*height; ++i) {
    if (pixel[i] == 0)
      pixel[i] = (int)(meanv-1.5*sdv);
  }

  // find the local standard deviation for each pixel
  int breakup = 10;
  float lmean;
  float smin = 10000.0;
  float smax = 0;
  int cc = 0;
  cnt = 0;

  for (int j = 0; j < width*height; ++j) {
    sddev[j] = 0.0;
    grid[j] = 0;
    grid2[j] = 0;
    locali[j] = 0;
  }

  // collect local standard deviation values for all pixels in the image
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
      int place = 0;
      for (int jj = firstj; jj < nextj; ++jj) {
	for (int ii = firstk; ii < nextk; ++ii) {
	  place = jj*width + ii;
	  lmean += (double)(pixel[place]);
	  ++cc;
	}
      }
      lmean = lmean/(double)(cc);
      locali[j*width+k] = (int)(lmean);
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

  cout << "sddev min and max: " << smin << " " << smax << endl;
  //cout << "current cnt: " << cnt << endl;

  // find the really bright pixels
  for (int j = 0; j < width*height; ++j) {
      bright2s[j] = 0;
      if (pixel[j] > (uint16)(meanv+2.0*sdv))
	bright2s[j] = 1;
  }

  // find locali clumps and store for colony separation later
  for (int j = 0; j < width*height; ++j) {
      grid[j] = 0;
      if (locali[j] > (int)(meanv+0.5*sdv))
	grid[j] = -1;
  }

  cc = floodfill(100);
  //showimage(grid); 
  //cout << "fcnt: " << cc << endl;
  //for (int j = 0; j < cc; ++j) {
    //cout << j << " " << clmpsize[j] << endl;
  //}

  // store the grid in locali
  for (int j = 0; j < width*height; ++j) {
    locali[j] = grid[j];
  }

  // use these clumps to divide colonies
  bool atedge;

  // the bright2s will show the shapes of the colonies: use the dark1s to
  // define clumps
  for (int j = 0; j < width*height; ++j) {
    grid[j] = 0;
    if (bright2s[j] > 0)
      grid[j] = 1;
  }
  dilate();
  dilate();
  dilate();

  for (int j = 0; j < width*height; ++j) {
    if (grid[j] > 0) grid[j] = -1;
  }
  cc = floodfill(100);
  //showimage(grid); 
  //cout << "fcnt: " << cc << endl;
  // store this grid in bright2s
  for (int j = 0; j < width*height; ++j) {
    bright2s[j] = grid[j];
  }

  // find roundness just to look at: store in grid2
  int rc = 0;
  float roundness;
  float circle;
  int extrap = 20;
  for (int j = 0; j < width*height; ++j) {
    grid2[j] = 0;
  }
  for (int p = 1; p <= cc; ++p) {
    if (clmpsize[p] > 2000) {
      //cout << "next clump: " << clmpsize[p] << " box: ";
      findthebox(p);
      //cout << bx1 << " " << bx2 << " " << by1 << " " << by2 << " ";
      rc = 0;
      for (int jj = by1-extrap; jj <= by2+extrap; ++jj) {
	for (int kk = bx1-extrap; kk <= bx2+extrap; ++kk) {
	  if (grid[jj*width+kk] > 0) ++rc;
	}
      }

      // measure roundness by circle that fits in box
      circle = pi*(float)(bx2-bx1+2*extrap)*(float)(bx2-bx1+2*extrap)/4.0;
      if ((by2-by1) > (bx2-bx1))
	circle = pi*(float)(by2-by1+2*extrap)*(float)(by2-by1+2*extrap)/4.0;
      roundness = (float)(rc)/circle;
      //cout << roundness << endl;
      // keep everything that is either big or round


      if (roundness < 0.2 || (clmpsize[p] < 500 && roundness < 0.3)) {      
	for (int jj = by1; jj <= by2; ++jj) {
	  for (int kk = bx1; kk <= bx2; ++kk) {
	    if (grid[jj*width+kk] == p) {
	      grid[jj*width+kk] = 0;
	      bright2s[jj*width+kk] = 0;
	    }
	  }
	}
      }

      // number in grid2 according to roundness      
      for (int jj = by1; jj <= by2; ++jj) {
	for (int kk = bx1; kk <= bx2; ++kk) {
	  if (grid[jj*width+kk] == p) {
	    grid2[jj*width+kk] = 1; //pink
	    if (roundness > 0.05) grid2[jj*width+kk] = 2; // red
	    if (roundness > 0.1) grid2[jj*width+kk] = 3; // orange
	    if (roundness > 0.2) grid2[jj*width+kk] = 4; // yellow
	    if (roundness > 0.3) grid2[jj*width+kk] = 5; // green
	    if (roundness > 0.4) grid2[jj*width+kk] = 6; // cyan
	    if (roundness > 0.5) grid2[jj*width+kk] = 7; // blue
	    if (roundness > 0.6) grid2[jj*width+kk] = 8; // purple
	    if (roundness > 0.7) grid2[jj*width+kk] = 9; // dark purple
	    if (roundness > 0.8) grid2[jj*width+kk] = 10; // light grey
	  }
	}
      }
    }
  }
  //showimage(grid2);

  // now dilate the bright pixels and fill holes
  for (int j = 0; j < width*height; ++j) {
    grid[j] = bright2s[j];
  }
  dilate();
  dilate();
  fillholes(5000); // fillholes uses grid and grid2

  // keep only the large clusters
  for (int j = 0; j < width*height; ++j) {
    if (grid[j] > 0)
      grid[j] = -1;
  }
  cc = floodfill(3000);
  //cout << "fcnt: " << cc << endl;
  //showimage(grid);

  int fr0 = 0;
  int fr1 = 0;
  int fr2 = 0;
  int fr3 = 0;
  int fr4 = 0;
  int fr5 = 0;
  int fr6 = 0;
  //float fractx[10];
  int bc = 0;
  float bavg = 0.0;
  for (int j = 0; j < width*height; ++j) {
    if (grid[j] == 0) {
      bavg += sddev[j];
      ++bc;
    }
  }
  bavg = bavg/(float)(bc);
  //printf("test against bavg: %f\n",bavg);
  bc = 0;
  float bsd = 1590.0;

  float sdmean = 0.0;
  float sdsd = 0.0;

  //float hg, contrast;
  colcnt = 0;
  int numclumps, lookfor;
  newcnt = cc;

  for (int j = 0; j < 10000; ++j) {
    perimx[j] = 0;
    perimy[j] = 0;
  }
  for (int j = 0; j < width*height; ++j) {
    perimeter[j] = 0;
  }

  for (int p = 1; p <= cc; ++p) {
    //printf(" %d %d\n",p,clmpsize[p]);

    findthebox(p);

    // process if close to center
    float cx = (float)(bx2+bx1)/2.0;
    float cy = (float)(by2+by1)/2.0;
    float dist = sqrt((cx-(width/2))*(cx-(width/2))+(cy-(height/2))*(cy-(height/2)));
    if (dist > 100.0) continue;
    atedge = false;
    if (bx1 < 10) atedge = true;
    if (bx2 > 1490) atedge = true;
    if (by1 < 10) atedge = true;
    if (by2 > 1410) atedge = true;
    if (atedge) continue;

    // for each, find how many pixels touch the edge
    bc = 0;
    fr1 = 0;
    fr2 = 0;
    fr3 = 0;
    fr4 = 0;
    fr5 = 0;
    fr6 = 0;
    ++colcnt;
    avgx = 0.0;
    avgy = 0.0;
    // process clump

    // is it only 1 clump: if 2, second clump returns with grid = -2, -3, etc.
    numclumps = isitone(p);
    //printf("numclumps from isitone: %d\n",numclumps);

    ///////////////////////////////////////////////////////////////////
    //numclumps = 1;
    ///////////////////////////////////////////////////////////////////
    //printf("box with possibly 2: %d %d %d %d\n",bx1,bx2,by1,by2);
    findthebox(p);
    //printf("box p: %d:    %d %d %d %d\n",p,bx1,bx2,by1,by2);
    findthebox(newcnt);
    //printf("new box p: %d:    %d %d %d %d\n",newcnt,bx1,bx2,by1,by2);
    for (int pp = 1; pp <= numclumps; ++pp) {
      lookfor = p;
      if (pp > 1) lookfor = newcnt;
      findthebox(lookfor);
      //printf("box: %d:    %d %d %d %d\n",pp,bx1,bx2,by1,by2);

      // smooth the edges, find the perimeter, fill
      if (pp == 1) smoothedges(lookfor);
      bc = 0;
      fr1 = 0;
      fr2 = 0;
      fr3 = 0;
      fr4 = 0;
      fr5 = 0;
      fr6 = 0;

      // remove border
      
      // count pixels in each sd range
      
      // define homogeneity and contrast
      int qq = 0;
      for (int jj = by1; jj <= by2; ++jj) {
	for (int kk = bx1; kk <= bx2; ++kk) {
	  if (grid[jj*width+kk] == lookfor) {
	    // add contribution to bsd calc if not a border pixels
	    if (grid[jj*width+kk-1] == lookfor && grid[jj*width+kk+1] == lookfor &&
		grid[(jj-1)*width+kk] == lookfor && grid[(jj+1)*width+kk] == lookfor &&
		grid[(jj-1)*width+kk-1] == lookfor && grid[(jj+1)*width+kk+1] == lookfor &&
		grid[(jj+1)*width+kk-1] == lookfor && grid[(jj-1)*width+kk+1] == lookfor &&
		grid[jj*width+kk-2] == lookfor && grid[jj*width+kk+2] == lookfor &&
		grid[(jj-2)*width+kk] == lookfor && grid[(jj+2)*width+kk] == lookfor &&
		grid[(jj-2)*width+kk-2] == lookfor && grid[(jj+2)*width+kk+2] == lookfor &&
		grid[(jj+2)*width+kk-2] == lookfor && grid[(jj-2)*width+kk+2] == lookfor &&
		grid[jj*width+kk-3] == lookfor && grid[jj*width+kk+3] == lookfor &&
		grid[(jj-3)*width+kk] == lookfor && grid[(jj+3)*width+kk] == lookfor &&
		grid[(jj-3)*width+kk-3] == lookfor && grid[(jj+3)*width+kk+3] == lookfor &&
		grid[(jj+3)*width+kk-3] == lookfor && grid[(jj-3)*width+kk+3] == lookfor) {
	      ++bc;
	      avgx += kk;
	      avgy += jj;
	      sdmean+= sddev[jj*width+kk];
	      if (sddev[jj*width+kk] <= bavg+0.5*bsd) ++fr0;
	      else if (sddev[jj*width+kk] < bavg+1.0*bsd) ++fr1;
	      else if (sddev[jj*width+kk] < bavg+1.5*bsd) ++fr2;
	      else if (sddev[jj*width+kk] < bavg+2.0*bsd) ++fr3;
	      else if (sddev[jj*width+kk] < bavg+2.5*bsd) ++fr4;
	      else if (sddev[jj*width+kk] < bavg+3.0*bsd) ++fr5;
	      else  ++fr6;
	      ++qq;
	    }
	  }
	}
      }
      //printf("count in main: %d\n",qq);
      if (bc > 0) {
	sdmean = sdmean/(float)(bc);      
	for (int jj = by1; jj <= by2; ++jj) {
	  for (int kk = bx1; kk <= bx2; ++kk) {
	    if (grid[jj*width+kk] == p) {
	      sdsd += (sddev[jj*width+kk]-sdmean)*(sddev[jj*width+kk]-sdmean);
	    }
	  }
	}
	sdsd = sqrt(sdsd/(float)(bc));
	avgx = avgx/(float)(bc);
	avgy = avgy/(float)(bc);
	centers[colcnt][0] = avgx;
	centers[colcnt][1] = avgy;
	//cout << "set centers: " << centers[colcnt][0] << " " << centers[colcnt][1] << endl;
	/*fractx[0] = (float)(fr0)/(float)(bc);
	fractx[1] = (float)(fr1)/(float)(bc);
	fractx[2] = (float)(fr2)/(float)(bc);
	fractx[3] = (float)(fr3)/(float)(bc);
	fractx[4] = (float)(fr4)/(float)(bc);
	fractx[5] = (float)(fr5)/(float)(bc);
	fractx[6] = (float)(fr6)/(float)(bc);
	contrast = fractx[6];
	hg = findhomogeneity(fractx);*/
	circularity = findtheperimeter(lookfor);
	//printf("circularity: %f\n",circularity);
	edgemeasure = findedgemeasure(lookfor);
      }
      
      // find the boundary of the colony and the center
      bx1 = bx1 - 50;
      if (bx1 < 0) bx1 = 0;
      bx2 = bx2 + 50;
      if (bx2 > width-1) bx2 = width - 1;
      by1 = by1 - 50;
      if (by1 < 0) by1 = 0;
      by2 = by2 + 50;
      if (by2 > height-1) by2 = height - 1;
      colboxes[colcnt][0] = bx1;
      colboxes[colcnt][1] = bx2;
      colboxes[colcnt][2] = by1;
      colboxes[colcnt][3] = by2;
      for (int jj = by1; jj <= by2; ++jj) {
	for (int kk = bx1; kk <= bx2; ++kk) {
	  perimeter[jj*width+kk] = 0;
	  if (grid[jj*width+kk] == p) {
	    if (grid[jj*width+kk-1] != p || grid[jj*width+kk+1] != p ||
		grid[(jj-1)*width+kk] != p || grid[(jj+1)*width+kk] != p) {
	      perimeter[jj*width+kk] = colcnt;
	    }
	  }
	}
      }
    }
    for (int h = 0; h < width*height; ++h) {
      if (perimeter[h] > 0)
	pixel[h] = 60000;
    }
    showimage2(pixel);
  }
}
