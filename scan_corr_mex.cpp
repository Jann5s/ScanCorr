#pragma warning( disable : 4267 4838 4244 4305) // cast warnings

/*
 [M, b, r, R*] = scan_corr_mex(f, g, a, bin, roi, mode*, interp*, Nthreads*)
 
 outputs (* are optional, can be empty)
 M    : sparse double ( 2N x 2N ) or ( 2N x 2 )
 b    : double ( 2N x 1 )
 r    : double ( 1 x 1 )
 R    : float ( n x m )
 
 inputs (* are optional, can be empty)
 f        : float ( n x m )
 g        : float ( n x m )
 a        : double ( 2N x 1 )
 bin      : double ( n x 1 ) or ( m x 1 ) example [ 0; 0; 0; 1; 1; 2; 2; 3; 3; ... n; n]
 roi      : double ( 2 x 1 ) [left, right] or [top, bottom] for mode 0, 1 resp.
 mode     : double ( 1 x 1 ) { 0, 1 } 0 = horizontal, 1 = vertical
 interp   : double ( 1 x 1 ) { 1 : 9 }
 Nthreads : double ( 1 x 1) { x >= 0, 0 = auto }
 
 */

// included for the mex code
#include "mex.h"

#include <string.h>
#include <limits>
#include <algorithm>

#include <thread>
#include <mutex>
#include <vector>

// include¡ for the bspline code
#include    "lib/interp_bspline.hpp"

// constants
enum Mode { horizontal, vertical };
enum Interp { i2 = 2, i3 = 3, i4 = 4, i5 = 5, i6 = 6, i7 = 7, i8 = 8, i9 = 9};
enum Flag { no, yes};
const char messageid[] = "scan_corr_mex:main";


// The stuct that contains a thread and all data
// ---------------------------------------------------------------------------
class ThreadData
{
    // std::mutex mutex;
public:

  // token for adding to global data
  std::mutex *mutex;
  int *visited_px;

  // list of pixels for the thread
  std::vector<int> list;
  
  // images
  float *f;
  float *g; // already as bsplineCoeff
  int n = 0;
  int m = 0;
  int N = 0;

	// region of interest
	int *roi; // [left, right, top, bottom]
  
  // interpolation
  int s = 0;
  int extent = 2;
  
  // dof and dic
  double *M;
  double *a;
  double *b;
  double *bin;

  double *r;
  float *R;
  
  // constructor
  ThreadData () {
  }
  
  // destructor
  ~ThreadData() {
  }
  
  // ----------------------------
  
  // some index operators
  int Ix(const int &I) const {
    return 2 * I;
  }
  int Iy(const int &I) const {
    return 2 * I + 1;
  }

  // image related functions
  // ----------------------------

  void ind2sub(int &i, int &j, const int &o) const {
    i = o % n;
    j = (o / n) % m;
  }

  void ind2sub(double &x, double &y, const int &o) const {
    int i, j;
    ind2sub(i, j, o);
    x = double(j);
    y = double(i);
  }

  int sub2ind(const int &i, const int &j) const {
    return i + n * j;
  }

	int is_inroi(const int &i) const {
		if ((i < roi[0]) || (i > roi[1] - 1)) {
			return 0;
		}
		else {
			return 1;
		}
	}

	int is_inroi(const double &x) const {
		return is_inroi(int(x));
	}

  int is_inside(const int &i, const int &j) const {
    if ((j < extent) || (j > m - 1 - extent) || (i < extent) || (i > n - 1 - extent)) {
			return 0;
		}
		else {
			return 1;
		}
  }
  
  int is_inside(const double &x, const double &y) const {
    if ((x < extent) || (x > m - 1 - extent) || (y < extent) || (y > n - 1 - extent)) {
			return 0;
		}
		else {
			return 1;
		}
	}

	 // interpolate on g
  float g_get_val(const double &x, const double &y) const {
    return (float)InterpolatedValue(g, n, m, y, x, s);
  }
  
  // get a value in f
  float f_get_val(const int &i, const int &j) const {
		int o = sub2ind(i, j);
		return f[o];
  }
  
  void f_get_grad(double *gra, const int &o) const {
    // at integer pixel location gradient
    int i, j;
    ind2sub(i, j, o);
    gra[0] = 0.5*(f_get_val(i + 0, j + 1) - f_get_val(i + 0, j - 1)); // (i, j) = (y, x) ==> dx
    gra[1] = 0.5*(f_get_val(i + 1, j + 0) - f_get_val(i - 1, j + 0)); // (i, j) = (y, x) ==> dy
    gra[0] = (isnan(gra[0])) ? 0 : gra[0];
    gra[1] = (isnan(gra[1])) ? 0 : gra[1];
  }

	void g_get_grad(double *gra, const double &x, const double &y) const {
		// at integer pixel location gradient
		gra[0] = 0.5*(g_get_val(x + 1, y + 0) - g_get_val(x - 1, y + 0)); // (i, j) = (y, x) ==> dx
		gra[1] = 0.5*(g_get_val(x + 0, y + 1) - g_get_val(x + 0, y - 1)); // (i, j) = (y, x) ==> dy
		// gra[0] = (isnan(gra[0])) ? 0 : gra[0];
		// gra[1] = (isnan(gra[1])) ? 0 : gra[1];
	}

  
  // DIC
  // ----------------------------
  template< Flag flagR , Mode mode >
  void dic(double *M_, double *b_, double &r_, int &cnt, const double *a, const int &o, const int &I) const {
    // compute and sum the local M, b, r for one pixel
    
		// get the reference coordinates
		int i, j;
		ind2sub(i, j, o);

		if (is_inside(i, j) == 0) {
			if (flagR == yes) {
				R[o] = std::numeric_limits<float>::quiet_NaN();
			}
			return;
		}

    // get the reference pixel
    float f_ = f[o];
      
      if (isnan(f_)){
          if (flagR == yes) {
              R[o] = std::numeric_limits<float>::quiet_NaN();
          }
          return;
      }

		// test the ROI
		if (mode == horizontal) {
			if (is_inroi(j) == 0) {
				if (flagR == yes) {
					R[o] = std::numeric_limits<float>::quiet_NaN();
				}
				return;
			}
		}
		else {
			if (is_inroi(i) == 0) {
				if (flagR == yes) {
					R[o] = std::numeric_limits<float>::quiet_NaN();
				}
				return;
			}
		}

		// get the deformed coordinates
		double x0 = int(j);
		double y0 = int(i);
		double x = x0 + a[Ix(I)];
    double y = y0 + a[Iy(I)];
    
    if (is_inside(x, y) == 0) {
        if (flagR == yes) {
            R[o] = std::numeric_limits<float>::quiet_NaN();
        }
        return;
    }
    
    // get the deformed pixel
    float g_ = g_get_val(x, y);
    
    // residual
    double res = double(f_) - double(g_);
    
    if (flagR == yes){
      R[o] = float(res);
    }
    
    // compute the gradient
    double grad[2];
    // g_get_grad(grad, x, y);
		f_get_grad(grad, o);
    
    // Hessian (only symmetric part)
    M_[0] += grad[0] * grad[0];
    M_[1] += grad[1] * grad[1];
    M_[2] += grad[0] * grad[1];

    // compute the right hand member
    b_[0] += grad[0] * res;
    b_[1] += grad[1] * res;

    // update the rms (only the squred sum, sqrt at the end)
    r_ += res*res;
    ++cnt;
  }
  
   
  // Thread execution
  // ----------------------------

  template< Flag flagR, Mode mode >
  std::thread spawn() {
    return std::thread([=] { exec<flagR, mode>(); } );
  }
  
  template< Flag flagR, Mode mode >
  void exec(){
    // create thread local variables

    // initiate the residual and counter for this thread
    double r_ = 0;
    int cnt = 0;
    
    // for each pixel row (in this thread)
    for (size_t l = 0; l < list.size(); l++) {
      
      // intiate the local M and b for this pixel row
      double M_[3] = { 0 };
      double b_[2] = { 0 };
      
      // start of the pixel row
      int o = list[l];
      const int I = int(bin[o]) - 1; // convert from ML style		
      
      // loop over all pixels in the row
      if (mode == horizontal) {
        for (int j = 0; j < m; j++){
          dic<flagR,mode>(M_, b_, r_, cnt, a, o + j * n, I);
        }
      }
      else if (mode == vertical) {
          for (int i = 0; i < n; i++){
          dic<flagR,mode>(M_, b_, r_, cnt, a, i + o * n, I);
        }
      }

      // lock
      mutex->lock();

      // write to global memory (this requires a lock if a bin is split over multiple threads)
      M[ Ix(I) + 2 * N * 0 ] += M_[0];
      M[ Iy(I) + 2 * N * 1 ] += M_[1];
      M[ Ix(I) + 2 * N * 1 ] += M_[2];
      M[ Iy(I) + 2 * N * 0 ] += M_[2];

      b[ Ix(I) ] += b_[0];
      b[ Iy(I) ] += b_[1];

      // unlock
      mutex->unlock();
    }
    
    // lock
    mutex->lock();

    r[0] += r_ ;
    visited_px[0] += cnt;

    // unlock
    mutex->unlock();
  }
};


// Main function (called by matlab)
// ---------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  // Defaults
  // ------------------------------------
  Interp s = static_cast<Interp>(3); // bspline mode
  double tol = 1e-9; // bspline tolerance
  Mode mode = horizontal; // 0 = H, 1 = V
  int Nthreads = 0; // 0 = auto

  // Processing intputs
  // ------------------------------------
  if (nrhs < 3) {
    mexErrMsgIdAndTxt(messageid,"Incorrect number of inputs inputs (too few).");
  }
  
	if (!mxIsClass(prhs[0],"single")) {
		mexErrMsgIdAndTxt(messageid, "gh must be single");
	}
  const int f_dims = mxGetNumberOfDimensions(prhs[0]);
  const mwSize *f_siz = mxGetDimensions(prhs[0]);
  float *f = (float *)mxGetData(prhs[0]) ;

	if (!mxIsClass(prhs[1], "single")) {
		mexErrMsgIdAndTxt(messageid, "gv must be single");
	}
	const int g_dims = mxGetNumberOfDimensions(prhs[1]);
  const mwSize *g_siz = mxGetDimensions(prhs[1]);
  float *g = (float *)mxGetData(prhs[1]) ;

	if (!mxIsClass(prhs[2], "double")) {
		mexErrMsgIdAndTxt(messageid, "a must be double");
	}
	const mwSize Na = mxGetNumberOfElements(prhs[2]);
  double *a = (double *)mxGetData(prhs[2]) ;

	if (!mxIsClass(prhs[3], "double")) {
		mexErrMsgIdAndTxt(messageid, "bin must be double");
	}
	const mwSize Nbin = mxGetNumberOfElements(prhs[3]);
  double *bin = (double *)mxGetData(prhs[3]);

	if (!mxIsClass(prhs[4], "double")) {
		mexErrMsgIdAndTxt(messageid, "roi must be double");
	}
	const mwSize Nroi = mxGetNumberOfElements(prhs[4]);
	double *roi_ = (double *)mxGetData(prhs[4]);

	if (Nroi != 2) {
		mexErrMsgIdAndTxt(messageid, "roi must be of size 2");
	}

  if ((nrhs > 5) && ( !mxIsEmpty( prhs[ 5 ] ) )){
		if (!mxIsClass(prhs[5], "double")) {
			mexErrMsgIdAndTxt(messageid, "mode must be double");
		}
		double mode_ = (double) mxGetScalar( prhs[ 5 ] );
    if (int(mode_) == 1) {
      mode = vertical;
    }
  }

  if ( ( nrhs > 6 ) && ( !mxIsEmpty( prhs[ 6 ] ) ) ) {
		if (!mxIsClass(prhs[6], "double")) {
			mexErrMsgIdAndTxt(messageid, "s must be double");
		}
		double s_ = (double) mxGetScalar( prhs[ 6 ] );
    if ( (int(s_) >= 2 ) && ( int(s_) <= 9) ) {
      s = static_cast<Interp>(int(s_));
    } else {
      mexErrMsgIdAndTxt(messageid,"The interp degree must be between 2 and 9 (inclusive).");
    }
  }

  if ( ( nrhs > 7 ) && ( !mxIsEmpty( prhs[ 7 ] ) ) ){
		if (!mxIsClass(prhs[7], "double")) {
			mexErrMsgIdAndTxt(messageid, "Nthreads must be double");
		}
		double Nthreads_ = (double) mxGetScalar( prhs[ 7 ] );
    if ( int(Nthreads_) >= 0 ) {
      Nthreads = int(Nthreads_);
    } else {
      mexErrMsgIdAndTxt(messageid,"The number of threads must be 0 or greater (0 = auto).");
    }
  }
  
  if (Nthreads == 0) {
      Nthreads = std::thread::hardware_concurrency();
  }

  // get some data sizes
  int n = f_siz[0];
  int m = f_siz[1];
  int N = int(Na/2);
  
  int roi[2] = { 0 };
	if (mode == horizontal) {
		roi[0] = std::max(int(roi_[0] - 1), 0);
		roi[1] = std::min(int(roi_[1] - 1), m);
	}
	else {
		roi[0] = std::max(int(roi_[0] - 1), 0);
		roi[1] = std::min(int(roi_[1] - 1), n);
	}
  
  for (int i = 1; i < 2 ; i++){
    if (f_siz[i] != g_siz[i]) {
      mexErrMsgIdAndTxt(messageid,"f and g must be the same size"); // TODO this is not required with a bit of smarter coding
    }
  }
  
  if ( (mode == horizontal) && (Nbin != n) ){
    mexErrMsgIdAndTxt(messageid,"The bin vector must match the image height (horizontal mode).");
  }
  if ( (mode == vertical) && (Nbin != m) ){
    mexErrMsgIdAndTxt(messageid,"The bin vector must match the image width (vertical mode).");
  }
  for (int i = 1; i < Nbin; i++){
    if ( bin[i] < 0) {
      mexErrMsgIdAndTxt(messageid,"bin values must be larger than zero");
    }
    if ( bin[i] - 1 >= N ) {
      mexErrMsgIdAndTxt(messageid,"bin values must be smaller than numel(a)/2");
    }
  }
    
  // Preparing outputs
  // ------------------------------------
  const Flag flagR = (nlhs == 4) ? yes : no;

  // M
  mwSize M_siz[2] = { mwSize(2*N), mwSize(2) };
  plhs[0] = mxCreateNumericArray(2,M_siz,mxDOUBLE_CLASS,mxREAL);
  double *M = ( (double *)mxGetData(plhs[0]) );

  // b
  mwSize b_siz[2] = { mwSize(2*N), mwSize(1) };
  plhs[1] = mxCreateNumericArray(2,b_siz,mxDOUBLE_CLASS,mxREAL);
  double *b = ( (double *)mxGetData(plhs[1]) );

  // r
  mwSize r_siz[2] = { mwSize(1), mwSize(1) };
  plhs[2] = mxCreateNumericArray(2,r_siz,mxDOUBLE_CLASS,mxREAL);
  double *r = ( (double *)mxGetData(plhs[2]) );

  // R
  float *R = nullptr;
  if (flagR == yes) {
    mwSize R_siz[2] = { mwSize(n), mwSize(m) };
    plhs[3] = mxCreateNumericArray(2, R_siz, mxSINGLE_CLASS, mxREAL);
    R = ((float *)mxGetData(plhs[3]));
  }

  // Preparing the Thread
  // ------------------------------------

  // convert g to Coeff data (for BSpline interp)
  float *c = new float[ n * m ];
  memcpy(c, g,  sizeof(float) * n * m );

	// make sure there are no nan's in the interpolated image
	for (int i = 0; i < n*m; i++) {
		if (isnan(c[i])) {
			c[i] = 0;
		}
	}

	// convert to BSpline coefficients
  int Error = SamplesToCoefficients(c, n, m, s, tol);

  // RMS counter
  int visited_px = 0;
  
  // create a mutex to control race conditions
  std::mutex mutex;
  
  // load the thread data
  ThreadData *TD = new ThreadData[ Nthreads ];

  for ( int i = 0; i < Nthreads ; i++) {
    TD[i].mutex = &mutex;
    
    // interpolation settings
    TD[i].s = static_cast<int>(s);
    TD[i].extent = std::max( ( TD[i].s + 2 ) / 2, 1);
    
    // image data
    TD[i].f = f;
    TD[i].g = c;
    TD[i].n = n;
    TD[i].m = m;

		TD[i].roi = roi;
    
    // DIC data
    TD[i].M = M;
    TD[i].a = a;
    TD[i].b = b;
    
    // shapefun
    TD[i].bin = bin;
    TD[i].N = N;
    
    TD[i].r = r;
    TD[i].R = R;
    TD[i].visited_px = &visited_px;
  }
  
  // count the number of rows that are active
  int Nrows = 0;
  for (int i = 0; i < Nbin ; i++) {
    if (bin[i] > 0) {
      Nrows++;
    }
  }

  // distribute the active pixel rows
  int cnt = 0;
  for (int o = 0; o < Nbin ; o++){
    if (bin[o] > 0){
      int t = int( double(Nthreads) * double(cnt) / double(Nrows) );
      TD[t].list.push_back( o );
      cnt++;
		}
		else {
			if (flagR == yes) {
				if (mode == horizontal) {
					for (int j = 0; j < m; j++) {
						R[o + j * n] = std::numeric_limits<float>::quiet_NaN();
					}
				}
				else if (mode == vertical) {
					for (int i = 0; i < n; i++) {
						R[i + o * n] = std::numeric_limits<float>::quiet_NaN();
					}
				}
			}
		}
  }
    
  // create a vector of threads
  std::vector<std::thread> v;
  
// the multithreaded work
	if (Nthreads > 1) {
		for (int i = 0; i < Nthreads; i++) {
			if ((flagR == yes) && (mode == horizontal)) {
				v.push_back(std::thread(TD[i].spawn<yes, horizontal>()));
			}
			else
				if ((flagR == yes) && (mode == vertical)) {
					v.push_back(std::thread(TD[i].spawn<yes, vertical>()));
				}
				else
					if ((flagR == no) && (mode == horizontal)) {
						v.push_back(std::thread(TD[i].spawn<no, horizontal>()));
					}
					else
						if ((flagR == no) && (mode == vertical)) {
							v.push_back(std::thread(TD[i].spawn<no, vertical>()));
						}
		}

		// combine the threads
		for (int i = 0; i < Nthreads; ++i) {
			v.at(i).join();
		}
	}
	else {
		for (int i = 0; i < Nthreads; i++) {
			if ((flagR == yes) && (mode == horizontal)) {
				TD[i].exec<yes, horizontal>();
			}
			else
				if ((flagR == yes) && (mode == vertical)) {
					TD[i].exec<yes, vertical>();
				}
				else
					if ((flagR == no) && (mode == horizontal)) {
						TD[i].exec<no, horizontal>();
					}
					else
						if ((flagR == no) && (mode == vertical)) {
							TD[i].exec<no, vertical>();
						}
		}
	}
	// finalize the rms
	r[0] = sqrt(r[0] / visited_px);

  // cleanup
  delete [] TD;
  delete [] c;

} // End of Main
