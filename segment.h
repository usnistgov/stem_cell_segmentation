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


class segment{
 public:
 segment(uint16* pix, int ww, int hh) {
   pixel = pix;
   width = ww;
   height = hh;
   sddev = new double[width*height];
   grid = new int[width*height];
   grid2 = new int[width*height];
   bright2s = new int[width*height];
   locali = new int[width*height];
   perimeter = new int[width*height];
   pedge2 = new int[width*height];
   qsize = 50000;
   maxcellsize = 20000;
 }

  ~segment();

  void findmask();
  int floodfill(int fmin);
  bool goodx(int x1);
  bool goody(int y1);
  void dilate();
  void dilateperim();
  void erodeperim();
  void findthebox(int p);
  void fillholes(int holesize);
  int isitone(int p);
  void smoothedges(int p);
  float findhomogeneity(float *fractx);
  float findtheperimeter(int lookfor);
  float findedgemeasure(int lookfor);
  void quickSort2(float* array, int* narray, int start, int end);
  void swap2(float* array, int* narray, int index1, int index2);

  void showimage(int* image);
  void showimage2(uint16* image);
  void showimage3(int* image);

  int width;
  int height;
  uint16* pixel;

  // need extra space to store pixel stddev values
  double *sddev;
  int *grid;
  int *grid2;
  int *locali;
  int* bright2s;
  int* perimeter;
  int* pedge2;
  int qsize;
  
  int clmpsize[50000];
  int queue[50000][2];
  int perimx[10000];
  int perimy[10000];
  float centers[50][2];
  int porder[1000];
  int colboxes[50][4];
  int steps[10000];
  int edgepix[100];
  int edgediffs[100000];
  float edgeratios[100000];

  float meanv, sdv;
  int bx1,bx2,by1,by2;
  int maxcellsize;
  int newcnt;
  int colcnt;
  float avgx, avgy;
  float circularity;
  float edgemeasure;
  int perimmark;
};
