#pragma warning( disable : 4267 4838 4244 4305) // cast warnings

/*
 f = scan_corr_mex(gh, gv, ah, ah, interp*, Nthreads*)
 
 outputs (* are optional, can be empty)
 f    : float ( n x m )
 
 inputs (* are optional, can be empty)
 gh       : float ( n x m )
 gv       : float ( n x m )
 ah       : double ( N x 2 )
 av       : double ( N x 2 )
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
enum Interp { i2 = 2, i3 = 3, i4 = 4, i5 = 5, i6 = 6, i7 = 7, i8 = 8, i9 = 9};
const char messageid[] = "scan_corr_im_mex:main";

// for dual threaded conversion
void BSCoeff(float *c, float *g, const int &n, const int &m, const int &s, const double &tol) {
	// convert gh and gv to Coeff data (for BSpline interp)
	memcpy(c, g, sizeof(float) * n * m);

	// make sure there are no nan's in the interpolated images
	for (int i = 0; i < n*m; i++) {
		if (isnan(c[i])) {
			c[i] = 0;
		}
	}

	// convert to BSpline coefficients
	int Error = SamplesToCoefficients(c, n, m, s, tol);
};


// The stuct that contains a thread and all data
// ---------------------------------------------------------------------------
class ThreadData
{
    // std::mutex mutex;
public:

  // list of pixels for the thread
  std::vector<int> list;
  
  // images
  float *gh; // already as bsplineCoeff
  float *gv; // already as bsplineCoeff
	float *f; // output image
	float *R; // residual image
	int n = 0;
  int m = 0;
  
  // interpolation
  int s = 0;
  int extent = 2;
  
  // dof
  double *ah;
	double *av;
  
  // constructor
  ThreadData () {
  }
  
  // destructor
  ~ThreadData() {
  }
  
  // ----------------------------

  // image related functions
  // ----------------------------

  void ind2sub(int &i, int &j, const int &o) const {
    i = o % n;
    j = (o / n) % m;
  }

  void ind2sub(double &x, double &y, const int &o) const {
    int i = o % n;
    int j = (o / n) % m;
    x = double(j);
    y = double(i);
  }

  int sub2ind(const int &i, const int &j) const {
    return i + n * j;
  }

  int is_inside(const int &x, const int &y) const {
    if ((x < extent) || (x > m - 1 - extent) || (y < extent) || (y > n - 1 - extent)) {
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

	 // interpolate on gh
  float gh_get_val(const float &x, const float &y) const {
    if (is_inside(x, y) == 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    return (float)InterpolatedValue(gh, n, m, y, x, s);
  }

	// interpolate on gv
	float gv_get_val(const float &x, const float &y) const {
		if (is_inside(x, y) == 0) {
			return std::numeric_limits<float>::quiet_NaN();
		}
		return (float)InterpolatedValue(gv, n, m, y, x, s);
	}

  // Merge
  // ----------------------------
  void merge(const int &o) const {
    // compute and sum the local M, b, r for one pixel
    
    // get the reference coordinates
    int i, j;
    ind2sub(i, j, o);

	  // get the H coordinates
    double xh = j + ah[i + 0 * n];
    double yh = i + ah[i + 1 * n];

		if (is_inside(xh, yh) == 0) {
			f[o] = std::numeric_limits<float>::quiet_NaN();
			return;
		}

		// get the V coordinates
		double xv = j + av[j + 0 * m];
		double yv = i + av[j + 1 * m];
    
		if (is_inside(xv, yv) == 0) {
			f[o] = std::numeric_limits<float>::quiet_NaN();
			return;
		}

		// get the deformed pixel
		float gh_ = gh_get_val(xh, yh);
		float gv_ = gv_get_val(xv, yv);

		// merge the images
		f[o] = 0.5 * (gh_ + gv_);

		// residual image
		R[o] = gh_ - gv_;
  }
  
  
  // Thread execution
  // ----------------------------

  std::thread spawn() {
    return std::thread([=] { exec(); } );
  }
  
	void exec() {
		// create thread local variables

		// for each pixel (in this thread)
		for (size_t l = 0; l < list.size(); l++) {

			// start of the pixel row
			int o = list[l];

			// merge the two images
			merge(o);
		}
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
  int Nthreads = 0; // 0 = auto

  // Processing intputs
  // ------------------------------------
  if (nrhs < 4) {
    mexErrMsgIdAndTxt(messageid,"Incorrect number of inputs inputs (too few).");
  }
  
	if (!mxIsClass(prhs[0], "single")) {
		mexErrMsgIdAndTxt(messageid, "gh must be single");
	}
	const int gh_dims = mxGetNumberOfDimensions(prhs[0]);
  const mwSize *gh_siz = mxGetDimensions(prhs[0]);
  float *gh = (float *)mxGetData(prhs[0]) ;

	if (!mxIsClass(prhs[1], "single")) {
		mexErrMsgIdAndTxt(messageid, "gh must be single");
	}
	const int gv_dims = mxGetNumberOfDimensions(prhs[1]);
  const mwSize *gv_siz = mxGetDimensions(prhs[1]);
  float *gv = (float *)mxGetData(prhs[1]) ;

	if (!mxIsClass(prhs[2], "double")) {
		mexErrMsgIdAndTxt(messageid, "ah must be double");
	}
	const int ah_dims = mxGetNumberOfDimensions(prhs[2]);
	const mwSize *ah_siz = mxGetDimensions(prhs[2]);
	double *ah = (double *)mxGetData(prhs[2]);
	
	if (!mxIsClass(prhs[3], "double")) {
		mexErrMsgIdAndTxt(messageid, "av must be double");
	}
	const int av_dims = mxGetNumberOfDimensions(prhs[3]);
	const mwSize *av_siz = mxGetDimensions(prhs[3]);
	double *av = (double *)mxGetData(prhs[3]);

  if ( ( nrhs > 4 ) && ( !mxIsEmpty( prhs[ 4 ] ) ) ) {
		if (!mxIsClass(prhs[4], "double")) {
			mexErrMsgIdAndTxt(messageid, "s must be double");
		}
		double s_ = (double) mxGetScalar( prhs[ 4 ] );
    if ( (int(s_) >= 2 ) && ( int(s_) <= 9) ) {
      s = static_cast<Interp>(int(s_));
    } else {
      mexErrMsgIdAndTxt(messageid,"The interp degree must be between 2 and 9 (inclusive).");
    }
  }

  if ( ( nrhs > 5 ) && ( !mxIsEmpty( prhs[ 5 ] ) ) ){
		if (!mxIsClass(prhs[5], "double")) {
			mexErrMsgIdAndTxt(messageid, "Nthreads must be double");
		}
		double Nthreads_ = (double) mxGetScalar( prhs[ 5 ] );
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
  int n = gh_siz[0];
  int m = gv_siz[1];

  for (int i = 1; i < 2 ; i++){
    if (gh_siz[i] != gv_siz[i]) {
      mexErrMsgIdAndTxt(messageid,"f and g must be the same size"); // TODO this is not required with a bit of smarter coding
    }
  }

	if (ah_siz[0] != n) {
		mexErrMsgIdAndTxt(messageid, "ah incorrect number of rows");
	}
	if (ah_siz[1] != 2) {
		mexErrMsgIdAndTxt(messageid, "ah incorrect number of cols");
	}
	if (av_siz[0] != m) {
		mexErrMsgIdAndTxt(messageid, "av incorrect number of rows");
	}
	if (av_siz[1] != 2) {
		mexErrMsgIdAndTxt(messageid, "av incorrect number of cols");
	}
    
  // Preparing outputs
  // ------------------------------------
  mwSize f_siz[2] = { mwSize(n), mwSize(m) };
  plhs[0] = mxCreateNumericArray(2, f_siz, mxSINGLE_CLASS, mxREAL);
  float *f = ((float *)mxGetData(plhs[0]));

	plhs[1] = mxCreateNumericArray(2, f_siz, mxSINGLE_CLASS, mxREAL);
	float *R = ((float *)mxGetData(plhs[1]));

  // Preparing the Thread
  // ------------------------------------

  // convert gh and gv to Coeff data (for BSpline interp)
  float *ch = new float[ n * m ];
	float *cv = new float[n * m];

	if (Nthreads > 1) {
		std::thread t1(BSCoeff, ch, gh, n, m, s, tol);
		std::thread t2(BSCoeff, cv, gv, n, m, s, tol);
		t1.join();
		t2.join();
	}
	else {
		BSCoeff(ch, gh, n, m, s, tol);
		BSCoeff(cv, gv, n, m, s, tol);
	}

  // load the thread data
  ThreadData *TD = new ThreadData[ Nthreads ];

  for ( int i = 0; i < Nthreads ; i++) {

    // interpolation settings
    TD[i].s = static_cast<int>(s);
		TD[i].extent = std::max((TD[i].s + 2) / 2, 1);
    
    // image data
    TD[i].f = f;
		TD[i].R = R;
		TD[i].gh = ch;
		TD[i].gv = cv;
		TD[i].n = n;
    TD[i].m = m;

		// dof
    TD[i].ah = ah;
    TD[i].av = av;    
  }

	// select a subset of pixel indices (simple segmentation, not taking the mask into account)
	for (int i = 0; i < Nthreads; i++) {
		int o1 = int(1.0 / double(Nthreads) * (i)* n*m);
		int o2 = int(1.0 / double(Nthreads) * (i + 1) * n*m);

		for (int o = o1; o < o2; o++) {
			TD[i].list.push_back(o);
		}
	}

  
    
  // create a vector of threads
  std::vector<std::thread> v;
  
  // the multithreaded work
  for ( int i = 0; i < Nthreads; i++) {
    v.push_back(std::thread(TD[i].spawn()));
  }
  
  // combine the threads
  for (int i=0; i < Nthreads; ++i){
    v.at(i).join();
  }

  // cleanup
  delete [] TD;
  delete [] ch;
	delete [] cv;

} // End of Main
