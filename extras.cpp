bool segment::goodx(int x1) {
  bool itsok = true;
  if (x1 < 0 || x1 > (width-1)) itsok = false;
  return itsok;
}

bool segment::goody(int y1) {
  bool itsok = true;
  if (y1 < 0 || y1 > (height-1)) itsok = false;
  return itsok;
}

int segment::floodfill(int fmin) {
  int cnt = 0;
  int clumpsize = 0;
  cnt = 0;
  int qsize = 50000;
  for (int j = 0; j < qsize; ++j) {
    queue[j][0] = 0;
    queue[j][1] = 0;
    clmpsize[j] = 0;
  }

  // if wholecells > 0, record box sizes
  int clpxmin, clpxmax, clpymin, clpymax;
  int bx1,bx2,by1,by2;
  for (int j = 1; j < height-1; ++j) {
    for (int k = 1; k < width-1; ++k) {
      // should you start here?
      //if (grid[j*width+k] == -1) cout << cnt << endl;
      if (grid[j*width+k] > -1) continue;
      ++cnt;
      bx1 = k;
      bx2 = k;
      by1 = j;
      by2 = j;
      //cout << "new point at:" << k << " " << j << endl;

      // new starting point
      int ptsqueue = 1;
      grid[j*width+k] = cnt;
      clumpsize = 1;
      queue[0][0] = k;
      queue[0][1] = j;
      clpxmin = k;
      clpxmax = k;
      clpymin = j;
      clpymax = j;
      while (ptsqueue > 0) {
	// find the next point to check
	int px = queue[0][0];
	int py = queue[0][1];

	// check surrounding
	if (goodx(px) && goody(py-1) && grid[(py-1)*width+px] == -1) {
	  //if (grid[px][py-1] == -1) {
	  //grid[px][py-1] = cnt;
	  grid[(py-1)*width+px] = cnt;
	  ++clumpsize;
	  if (px < clpxmin) clpxmin = px;
	  if (px > clpxmax) clpxmax = px;
	  if ((py-1) < clpymin) clpymin = py-1;
	  if ((py-1) > clpymax) clpymax = py-1;
	  // add point to queue
	  queue[ptsqueue][0] = px;
	  queue[ptsqueue][1] = py-1;
	  ++ptsqueue;
	  if (px < bx1) bx1 = px;
	  if (px > bx2) bx2 = px;
	  if ((py-1) < by1) by1 = py - 1;
	  if ((py-1) > by2) by2 = py - 1;
	}
	if (goodx(px+1) && goody(py) && grid[py*width+(px+1)] == -1) {
	  //if (grid[px+1][py] == -1) {
	  //grid[px+1][py] = cnt;
	  grid[py*width+(px+1)] = cnt;
	  ++clumpsize;
	  if ((px+1) < clpxmin) clpxmin = px+1;
	  if ((px+1) > clpxmax) clpxmax = px+1;
	  if (py < clpymin) clpymin = py;
	  if (py > clpymax) clpymax = py;
	  // add point to queue
	  queue[ptsqueue][0] = px+1;
	  queue[ptsqueue][1] = py;
	  ++ptsqueue;
	  if ((px+1) < bx1) bx1 = px + 1;
	  if ((px+1) > bx2) bx2 = px + 1;
	  if (py < by1) by1 = py;
	  if (py > by2) by2 = py;
	}
	if (goodx(px) && goody(py+1) && grid[(py+1)*width+px] == -1) {
	  //if (grid[px][py+1] == -1) {
	  //grid[px][py+1] = cnt;
	  grid[(py+1)*width+px] = cnt;
	  ++clumpsize;
	  if (px < clpxmin) clpxmin = px;
	  if (px > clpxmax) clpxmax = px;
	  if ((py+1) < clpymin) clpymin = py+1;
	  if ((py+1) > clpymax) clpymax = py+1;
	  // add point to queue
	  queue[ptsqueue][0] = px;
	  queue[ptsqueue][1] = py+1;
	  ++ptsqueue;
	  if (px < bx1) bx1 = px;
	  if (px > bx2) bx2 = px;
	  if ((py+1) < by1) by1 = py + 1;
	  if ((py+1) > by2) by2 = py + 1;
	}
	if (goodx(px-1) && goody(py) && grid[py*width+(px-1)] == -1) {
	  //if (grid[px-1][py] == -1) {
	  //grid[px-1][py] = cnt;
	  grid[py*width+(px-1)] = cnt;
	  ++clumpsize;
	  if ((px-1) < clpxmin) clpxmin = px-1;
	  if ((px-1) > clpxmax) clpxmax = px-1;
	  if (py < clpymin) clpymin = py;
	  if (py > clpymax) clpymax = py;
	  // add point to queue
	  queue[ptsqueue][0] = px-1;
	  queue[ptsqueue][1] = py;
	  ++ptsqueue;
	  if ((px-1) < bx1) bx1 = px - 1;
	  if ((px-1) > bx2) bx2 = px - 1;
	  if (py < by1) by1 = py;
	  if (py > by2) by2 = py;
	}

	// take this point off the queue
	for (int ii = 1; ii < ptsqueue; ++ii) {
	  queue[ii-1][0] = queue[ii][0];
	  queue[ii-1][1] = queue[ii][1];
	}
	--ptsqueue;
      }

      // check the clumpsize
      //cout << "clumpsize: " << clumpsize << endl;

      // if too small a clump: zero out the grid
      if (clumpsize < fmin) {
	for (int jj = clpymin; jj <= clpymax; ++jj) {
	  for (int kk = clpxmin; kk <= clpxmax; ++kk) {
	    if (grid[jj*width+kk] == cnt)
	      grid[jj*width+kk] = 0;
	  }
	}
	cnt = cnt - 1;
      }
      else {
	//printf("next clump: %d\n",clumpsize);
	clmpsize[cnt] = clumpsize;
      }
    }
  }
  //printf("cnt2: %d\n",cnt);
  int finalcnt = cnt;
  return (finalcnt);
}

void segment::smoothedges(int p) {
  //cout << "start of smoothdeges: " << p << endl;
  // model each colony as a circle or ellipse and fill in dark edges
  float ratio = (float)(by2-by1)/(float)(bx2-bx1);
  //printf("next ratio: %f\n",ratio);
  //int slices = 100;
  //float eachslice = 2.0*pi/(float)(slices);
  //printf("each slice: %f\n",eachslice);
  float dist, npx, npy;
  float cx = (float)(bx2+bx1)/2.0;
  float cy = (float)(by2+by1)/2.0;
  float radius = (float)(bx2-bx1)/2.0;
  if ((by2-by1) < (bx2-bx1))
    radius = (float)(by2-by1)/2.0;
  //printf("radius: %f\n",radius);
  //if (ratio > 0.75 && ratio < 1.75) {
  if (ratio > 0.65 && ratio < 1.75) {
    // circle
    //printf("circle\n");
    // find the perimeter
    for (int jj = by1; jj <= by2; ++jj) {
      for (int kk = bx1; kk <= bx2; ++kk) {
	perimeter[jj*width+kk] = 0;
	if (grid[jj*width+kk] == p) {
	  npx = (float)(kk);
	  npy = (float)(jj);
	  dist = sqrt((cx-npx)*(cx-npx)+(cy-npy)*(cy-npy));
	  if (fabs(dist-radius) < radius/2.0) {
	    perimeter[jj*width+kk] = (int)(dist);
	  }
	}
      }
    }
  }
  else {
    // ellipse
    //printf("ellipse\n");
  }
  //showimage3(perimeter);
}

void segment::erodeperim() {
  // erode the pattern in grid, copy into testi
  for (int j = 0; j < width*height; ++j) {
    grid2[j] = 0;
  }

  for (int j = 0; j < height; ++j) {
    for (int k = 0; k < width; ++k) {
      if (perimeter[j*width+k] > 0) {
	// keep if everything around it is also > 0
	if (perimeter[(j-1)*width+k-1] > 0 &&
	    perimeter[j*width+k-1] > 0 &&
	    perimeter[(j+1)*width+k-1] > 0 &&
	    perimeter[(j-1)*width+k] > 0 &&
	    perimeter[(j+1)*width+k] > 0 &&
	    perimeter[(j-1)*width+k+1] > 0 &&
	    perimeter[j*width+k+1] > 0 &&
	    perimeter[(j+1)*width+k+1] > 0)
	  grid2[j*width+k] = 1;
      }
    }
  }

  for (int j = 0; j < width*height; ++j) {
      perimeter[j] = grid2[j];
  }
}

void segment::dilate() {
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

void segment::dilateperim() {
  // dilate the pattern in grid, copy into testi

  //printf("box in perimdilate: %d %d %d %d\n",bx1,bx2,by1,by2);
  //int cnt = 0;
  for (int j = by1-5; j <= by2+5; ++j) {
    for (int k = bx1-5; k <= bx2+5; ++k) {
      grid2[j*width+k] = perimeter[j*width+k];
      //if (perimeter[j*width+k] > 0) ++cnt;
    }
  }
  //cout << "cnt: " << cnt << endl;

  for (int j = by1; j <= by2; ++j) {
    for (int k = bx1; k <= bx2; ++k) {
      if (perimeter[j*width+k] > 0) {
	grid2[j*width+k] = 1;
	grid2[j*width+k-1] = 1;
	grid2[j*width+k+1] = 1;
	grid2[(j-1)*width+k] = 1;
	grid2[(j+1)*width+k] = 1;
      }
    }
  }

  //cnt = 0;
  for (int j = by1; j <= by2; ++j) {
    for (int k = bx1; k <= bx2; ++k) {
      if (grid2[j*width+k] > 0) {
	perimeter[j*width+k] = 100;
	//++cnt;
      }
    }
  }
  //cout << "cnt2: " << cnt << endl;
}

void segment::findthebox(int p) {
  // find the box that contains the whole cell
  bx1 = 10000;
  bx2 = -1;
  by1 = 10000;
  by2 = -1;
  for (int j = 0; j < height; ++j) {
    for (int k = 0; k < width; ++k) {
      if (grid[j*width+k] == p) {
	if (k < bx1) bx1 = k;
	if (k > bx2) bx2 = k;
	if (j < by1) by1 = j;
	if (j > by2) by2 = j;
      }
    }
  }
  //printf("box for %d: size: %d  %d %d %d %d\n",p,clmpsize[p],bx1,bx2,by1,by2);
}

void segment::fillholes(int holesize) {
  // fill small holes in grid array

  // first store in grid2
  for (int j = 0; j < width*height; ++j) {
    grid2[j] = grid[j];
  }

  for (int j = 0; j < height; ++j) {
    for (int k = 0; k < width; ++k) {
      if (grid[j*width+k] == 0)
	 grid[j*width+k] = -1;
      else
	grid[j*width+k] = 0;
      // fill in the boundaries
      if (k == 0 || k == (width-1) || j == 0 || j == (height-1))
	grid[j*width+k] = 0;
    }
  }
  int cc = floodfill(1);
  //printf("called floodfill: %d\n",cc);

  for (int p = 1; p < cc; ++p) {
    if (clmpsize[p] < holesize) {
      for (int j = 0; j < width*height; ++j) {
	if (grid[j] == p)
	  grid2[j] = 1;
      }
    }
  }

  for (int j = 0; j < width*height; ++j) {
    grid[j] = grid2[j];
  }
}

int segment::isitone(int p) {
  int newnumber = 1;
  // look at grid = p to divide the colony
  int st, ed, lowpt, lowsize;
  int skipit;
  float fr;
  bool issplit;
  int cnt1, cnt2;
  
  // see if the clump is circular
  skipit = (by2 - by1)/8;
  for (int jj = by1; jj <= by2; ++jj) {
    perimy[jj] = 0;
    // look at each row: store (last-first) in perimy
    st = bx2;
    ed = bx1;
    for (int kk = bx1; kk <= bx2; ++kk) {
      if (grid[jj*width+kk] == p) {
	if (kk < st) st = kk;
	if (kk > ed) ed = kk;
      }
    }
    perimy[jj] = ed - st + 1;
    //if (p == 4) printf("next jj: %d diff: %d\n",jj,perimy[jj]);
  }

  // find lowest
  lowpt = by2;
  lowsize = by2 - by1;
  for (int jj = by1+skipit; jj <= by2-skipit; ++jj) {
    if (perimy[jj] < lowsize) {
      lowsize = perimy[jj];
      lowpt = jj;
    }
  }
  // compute what fraction the low pt is at
  fr = (float)(lowpt-by1-skipit)/(float)(by2-by1-2*skipit);
  //printf("lowpty: %d lowsize: %d  fraction: %f\n",lowpt,lowsize,fr);
  issplit = false;
  if (fr > 0.25 && fr < 0.75) {
    cnt1 = 0;
    cnt2 = 0;
    for (int jj = by1; jj <= by2; ++jj) {
      for (int kk = bx1; kk <= bx2; ++kk) {
	if (grid[jj*width+kk] == p) {
	  if (jj < lowpt) ++cnt1;
	  else ++cnt2;
	}
      }
    }
    //printf("each piece split in x direction: %d %d\n",cnt1,cnt2);
    
    // split this colony into 2 if both pieces are big enough
    if (cnt1 < maxcellsize) {
      //printf("first half x too small\n");
      for (int jj = by1; jj <= by2; ++jj) {
	for (int kk = bx1; kk <= bx2; ++kk) {
	  if (grid[jj*width+kk] == p) {
	    if (jj < lowpt)
	      grid[jj*width+kk] = 0;
	  }
	}
      }
    }
    if (cnt2 < maxcellsize) {
      //printf("second half x too small\n");
      for (int jj = by1; jj <= by2; ++jj) {
	for (int kk = bx1; kk <= bx2; ++kk) {
	  if (grid[jj*width+kk] == p) {
	    if (jj >= lowpt)
	      grid[jj*width+kk] = 0;
	  }
	}
      }
    }
    if (cnt1 > maxcellsize && cnt2 > maxcellsize) {
      newcnt = newcnt + 1;
      newnumber = newnumber + 1;
      //printf("increase newcnt a: %d\n",newcnt);
      issplit = true;
      for (int jj = by1; jj <= by2; ++jj) {
	for (int kk = bx1; kk <= bx2; ++kk) {
	  if (grid[jj*width+kk] == p) {
	    if (jj < lowpt) {
	    }
	    else {
	      grid[jj*width+kk] = newcnt;
	    }
	  }
	}
      }
    }
  }
  
  if (!issplit) {
    skipit = (bx2 - bx1)/8;
    for (int kk = bx1; kk <= bx2; ++kk) {
      perimx[kk] = 0;
      // look at each column: store (last-first) in perimx
      st = by2;
      ed = by1;
      for (int jj = by1; jj <= by2; ++jj) {
	if (grid[jj*width+kk] == p) {
	  if (jj < st) st = jj;
	  if (jj > ed) ed = jj;
	}
      }
      perimx[kk] = ed - st + 1;
      //printf("next kk: %d diff: %d\n",kk,perimx[kk]);
    }
    
    // find lowest
    lowpt = bx2;
    lowsize = bx2 - bx1;
    for (int kk = bx1+skipit; kk <= bx2-skipit; ++kk) {
      if (perimx[kk] < lowsize) {
	lowsize = perimx[kk];
	lowpt = kk;
      }
    }
    fr = (float)(lowpt-bx1-skipit)/(float)(bx2-bx1-2*skipit);
    //printf("lowptx: %d lowsize: %d  fraction: %f\n",lowpt,lowsize,fr);
    if (fr > 0.25 && fr < 0.75) {
      cnt1 = 0;
      cnt2 = 0;
      for (int jj = by1; jj <= by2; ++jj) {
	for (int kk = bx1; kk <= bx2; ++kk) {
	  if (grid[jj*width+kk] == p) {
	    if (kk < lowpt) ++cnt1;
	    else ++cnt2;
	  }
	}
      }
      //printf("each piece: %d %d\n",cnt1,cnt2);
      
      // split this colony into 2 if both pieces are big enough
      if (cnt1 < maxcellsize) {
	//printf("first half y too small\n");
	for (int jj = by1; jj <= by2; ++jj) {
	  for (int kk = bx1; kk <= bx2; ++kk) {
	    if (grid[jj*width+kk] == p) {
	      if (kk < lowpt)
		grid[jj*width+kk] = 0;
	    }
	  }
	}
      }
      if (cnt2 < maxcellsize) {
	//printf("second half y too small\n");
	for (int jj = by1; jj <= by2; ++jj) {
	  for (int kk = bx1; kk <= bx2; ++kk) {
	    if (grid[jj*width+kk] == p) {
	      if (kk >= lowpt)
		grid[jj*width+kk] = 0;
	    }
	  }
	}
      }
      if (cnt1 > maxcellsize && cnt2 > maxcellsize) {
	// split this colony into 2
	newcnt = newcnt + 1;
	newnumber = newnumber + 1;
	//printf("increase newcnt b: %d\n",newcnt);
	for (int jj = by1; jj <= by2; ++jj) {
	  for (int kk = bx1; kk <= bx2; ++kk) {
	    if (grid[jj*width+kk] == p) {
	      if (kk < lowpt) {
	      }
	      else {
		grid[jj*width+kk] = newcnt;
	      }
	    }
	  }
	}
      }
    }
  }

  return(newnumber);
}

float segment::findhomogeneity(float *fractx) {
  // order them, then see how many to get to 0.75
  fractx[5] += fractx[6];
  for (int pp = 0; pp < 6; ++pp) {
    // for (int pp = 0; pp < 7; ++pp) {
    porder[pp] = pp;
  }
  //printf("before: %f %f %f %f %f %f %f\n",fractx[0],fractx[1],fractx[2],fractx[3],fractx[4],fractx[5],fractx[6]);
  quickSort2(fractx, porder, 0, 5);
  //printf("after: %f %f %f %f %f %f %f\n",fractx[0],fractx[1],fractx[2],fractx[3],fractx[4],fractx[5],fractx[6]);

  float sum = fractx[5]+0.5*fractx[4];

  // min value if all equal = 1/6 * 3/2 = .25: so subtract .25 and divide by .75
  sum = (sum - 0.25)/0.75;
  return(sum);
}


float segment::findtheperimeter(int lookfor) {
  int pc = 0;
  int gridcnt = 0;
  //printf("limits: %d %d %d %d\n",bx1,bx2,by1,by2);
  for (int jj = by1; jj <= by2; ++jj) {
    for (int kk = bx1; kk <= bx2; ++kk) {
      perimeter[jj*width+kk] = 0;
      if (grid[jj*width+kk] == lookfor) {
	++gridcnt;
	if (grid[jj*width+kk-1] == 0 || grid[jj*width+kk+1] == 0 ||
	    grid[(jj-1)*width+kk] == 0 || grid[(jj+1)*width+kk] == 0) {
	  perimeter[jj*width+kk] = 1;
	  ++pc;
	}
      }
    }
  }
  //printf("start of findtheperimeter: pc: %d gridcnt: %d\n",pc,gridcnt);
  dilateperim();

  // try something new: split into 1000 angles
  float angle;
  float dist, dx, dy, posx, posy;
  int nsteps;
  int numpoints = 4000;
  for (int i = 0; i < numpoints; ++i) {
    perimx[i] = 0;
    perimy[i] = 0;
    steps[i] = 0;
  }
  int perimcnt = 0;
  //printf("start at center: %f %f\n",avgx,avgy);
  for (int ang = 0; ang < numpoints; ++ang) {
    angle = 2 * pi * (float)(ang) / (float)(numpoints);
    dx = cos(angle);
    dy = sin(angle);
    posx = avgx;
    posy = avgy;
    nsteps = 0;
    while (nsteps < 500) {
      posx = posx + 2.0*dx;
      posy = posy + 2.0*dy;
      if (perimeter[(int)(posy)*width+(int)(posx)] > 0) {
	perimx[perimcnt] = (int)(posx);
	perimy[perimcnt] = (int)(posy);
	steps[perimcnt] = nsteps;
	//printf("next angle: %f  cos: %f sin: %f pos: %f %f nsteps: %d\n",angle,cos(angle),sin(angle),posx,posy,nsteps);
	++perimcnt;
	break;
      }
      ++nsteps;
    }
  }

  circularity = 0.0;

  float sum = 0.0;
  int usept[numpoints];
  for (int w = 0; w < numpoints; ++w) {
    usept[w] = 0;
  }

  // use every 100th point
  int pickit = 50;
  int minsteps = 10000;
  int useit;
  int next = 0;
  for (int w = 0; w < perimcnt/pickit; ++w) {
    minsteps = 10000;
    useit = pickit*w;
    for (int ww = w*pickit; ww < w*pickit+pickit; ++ww) {
      if (steps[ww] < minsteps) {
	minsteps = steps[ww];
	useit = ww;
      }
    }
    //printf("use this one: useit: %d minsteps: %d\n",useit,minsteps);
    usept[next] = useit;
    ++next;
  }
  //printf("next at end: %d\n",next);

  float meanr = 0.0;
  for (int w = 1; w < next; ++w) {
    dx = perimx[usept[w]] - perimx[usept[w-1]];
    dy = perimy[usept[w]] - perimy[usept[w-1]];
    dist = sqrt(dx*dx + dy*dy);
    sum += dist;
    meanr += (float)(steps[usept[w]]);
  }
  meanr = meanr/(float)(next);
  float sdr = 0.0;
  for (int w = 1; w < next; ++w) {
    sdr += (meanr-(float)(steps[usept[w]]))*(meanr-(float)(steps[usept[w]]));
  }
  sdr = sqrt(sdr/(float)(next));
  //printf("meanr, sdr: %f %f perimcnt2: %d\n",meanr,sdr,perimcnt);

  // add first to last
  dx = perimx[usept[0]] - perimx[usept[next-1]];
  dy = perimy[usept[0]] - perimy[usept[next-1]];
  dist = sqrt(dx*dx + dy*dy);
  sum += dist;
  //printf("new count: %d final sum: %f\n",perimcnt,sum);

  perimmark = 100;
  for (int jj = by1; jj <= by2; ++jj) {
    for (int kk = bx1; kk <= bx2; ++kk) {
      perimeter[jj*width+kk] = 0;
      pedge2[jj*width+kk] = 0;
    }
  }

  for (int w = 1; w < perimcnt; ++w) {
    int kk = perimx[w];
    int jj = perimy[w];
    perimeter[jj*width+kk] = 1;
    pedge2[jj*width+kk] = 1;
  }

  // dilate perimeter for edge calculation: find pixels inside perimeter
  dilateperim();
  /*for (int jj = by1; jj <= by2; ++jj) {
    for (int kk = bx1; kk <= bx2; ++kk) {
      if (grid[jj*width+kk] == lookfor) {
	// zero out if 3 from edge
	grid[jj*width+kk] = -1;
	if (perimeter[jj*width+kk] > 0)
	  grid[jj*width+kk] = 0;
      }
    }
    }*/
  /*int cc = floodfillpart(1000);
  // pull out biggest
  for (int p = 1; p <= cc; ++p) {
    printf("pieces: %d\n",clmpsize2[p]);
    findtheboxgrid3(p);
    if (cx1 > 1 && cy1 > 1)
      bc = clmpsize2[p];
      }*/

  circularity = 4.0 * pi * (float)(gridcnt);
  circularity = circularity/(sum*sum);
  //printf("circularity: %f area: %d sum: %f\n",circularity,gridcnt,sum);

  // findegemeasure is expecting a dilated perimeter
  for (int times = 0; times < 20; ++times) {
    dilateperim();
  }

  int newcnt = 0;
  for (int w = 1; w < perimcnt; ++w) {
    int kk = perimx[w];
    int jj = perimy[w];
    perimeter[jj*width+kk] = perimmark;
    pedge2[jj*width+kk] = perimmark;
    ++newcnt;
  }
  //printf("end of circularity newcnt: %d\n",newcnt);
  return(circularity);
}

float segment::findedgemeasure(int lookfor) {
  // the inside of the perimeter rim is marked
  // for each marked perimeter pixel, find intensities outward
  int cc = 0;
  int rimcnt;
  float dx, dy, posx, posy;
  int minpix, maxpix;
  for (int jj = by1; jj <= by2; ++jj) {
    for (int kk = bx1; kk <= bx2; ++kk) {
      if (perimeter[jj*width+kk] == perimmark) {
	rimcnt = 0;
	posx = (float)(kk);
	posy = (float)(jj);
	dx = posx - avgx;
	dy = posy - avgy;
	dx = dx/sqrt(dx*dx + dy*dy);
	dy = dy/sqrt(dx*dx + dy*dy);
	//if (cc < 3) printf("center of colony: %f %f   kk,jj: %d %d  dx,dy: %f %f\n",avgx,avgy,kk,jj,dx,dy);
	// step further out until distance > cb or you hit a pixel with perimeter = 0
	while (rimcnt < 50) {
	  posx = posx + 2.0*dx;
	  posy = posy + 2.0*dy;
	  //if (cc < 3) printf("next point: %f %f\n",posx,posy);
	  if (perimeter[(int)(posy)*width+(int)(posx)] == 0)
	    break;
	  else {
	    edgepix[rimcnt] = pixel[(int)(posy)*width+(int)(posx)];
	    ++rimcnt;
	  }
	}
	if (rimcnt > 0) {
	  minpix = 100000;
	  maxpix = 0;
	  // find biggest gradient
	  for (int w = 0; w < rimcnt; ++w) {
	    if (edgepix[w] < minpix) minpix = edgepix[w];
	    if (edgepix[w] > maxpix) maxpix = edgepix[w];
	  }
	  //if (cc < 10) printf("rimcnt: %d min and maxpix; %d %d\n",rimcnt,minpix,maxpix);
	  edgediffs[cc] = maxpix - minpix;

	  // this time make it ratio of outside avg/inside avg
	  if (rimcnt/2 > 1) {
	    minpix = 100000;
	    maxpix = 0;
	    for (int w = 0; w < rimcnt/2; ++w) {
	      if (edgepix[w] > maxpix) maxpix = edgepix[w];
	      if (edgepix[rimcnt-1-w] < minpix) minpix = edgepix[rimcnt-1-w];
	    }
	    //edgeratios[cc] = (float)(maxpix)/(float)(minpix);
	    // make edge measurement a function of stddev's apart
	    edgeratios[cc] = (float)(maxpix-minpix)/sdv;
	    //if (cc < 10) printf("rimcnt: %d min and max; %d %d  ratio: %f\n",rimcnt,minpix,maxpix,edgeratios[cc]);
	  }
	  else
	    --cc;
	}
	++cc;
      }
    }
  }
  printf("final cc: %d final rimcnt: %d\n",cc,rimcnt);

  // find average
  float avgp = 0.0;
  for (int w = 0; w < cc; ++w) {
    avgp += edgeratios[w];
    //printf("edgeratio: %f avgp: %f\n",edgeratios[w],avgp);
  }
  edgemeasure = avgp/(float)(cc);
  //printf("edgemeasure: %f\n",edgemeasure);
  // 5 is the most
  edgemeasure = edgemeasure/5.0;
  return (edgemeasure);
}


void segment::swap2(float* array, int* narray, int index1, int index2) {
  float temp = array[index1];
  array[index1] = array[index2]; 
  array[index2] = temp;
  int tempi = narray[index1];
  narray[index1] = narray[index2]; 
  narray[index2] = tempi;
}

void segment::quickSort2(float* array, int* narray, int start, int end) {
  int i = start;
  int k = end;

  if (end - start >= 1)  {
    float pivot = array[start]; 

    while (k > i) {
      while (array[i] <= pivot && i <= end && k > i)
	i++; 
      while (array[k] > pivot && k >= start && k >= i)
	k--; 
      if (k > i)
	swap2(array, narray, i, k);
    }
    swap2(array, narray, start, k);
    quickSort2(array, narray, start, k - 1);
    quickSort2(array, narray, k+1, end);
  }
  else  {
    return;
  }
}

void segment::showimage2(uint16* image) {
  TIFF *tif = TIFFOpen("/home/peskin/compbio/leiber/wipp/image.tif","w");
  uint16 *buf = new uint16[width];
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; ++j) {
      buf[j] = (uint16)(image[width*i + j]);
    }
    TIFFWriteScanline(tif, buf, i, 0);
  }
  TIFFClose(tif);
  delete [] buf;
}

void segment::showimage(int* image) {
  TIFF *tif = TIFFOpen("/home/peskin/compbio/leiber/wipp/image.tif","w");
  uint8 *buf = new uint8[width];
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; ++j) {
      buf[j] = (uint8)(image[width*i + j]);
    }
    TIFFWriteScanline(tif, buf, i, 0);
  }
  TIFFClose(tif);
  delete [] buf;
}

void segment::showimage3(int* image) {
  TIFF *tif = TIFFOpen("/home/peskin/compbio/leiber/wipp/image3.tif","w");
  uint8 *buf = new uint8[width];
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; ++j) {
      buf[j] = (uint8)(image[width*i + j]);
    }
    TIFFWriteScanline(tif, buf, i, 0);
  }
  TIFFClose(tif);
  delete [] buf;
}
