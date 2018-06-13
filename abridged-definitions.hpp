#ifndef GAIA_DEFINITIONS
#define GAIA_DEFINITIONS
#define GAIA_VERSION 0.903

//Un-comment next #define to compile with GSL. Comment if GSL ISN'T installed
//and change LFLAGS in Makefile so it doesn't try to link ANY GSL libs.
//#define GSL_HERE

//Un-comment next #define to compile with NetCDF. Comment if NetCDF ISN'T installed
//and change LFLAGS in Makefile so it doesn't try to link ANY NetCDF libs.
//#define NETCDF_HERE

//Un-comment next #define to compile with NFFT, NFFTLS. Comment if NFFT ISN'T installed
//and change LFLAGS in Makefile so it doesn't try to link ANY NFFT libs.
//#define NFFT_HERE

//Un-comment next #define if GMT v5 is used (this isn't finished). Comment if GMT v4 is used.
//#define GMT5

//Load C libraries which are common to nearly all GAIA programs:
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>//For realpath().
//Allows the use of cos() and sin(), among other functions:
#include <math.h>
//Used for searching directories:
#include <dirent.h>
//Used for sizing notification intervals in loops on the fly and
//for providing estimates of time remaining:
#include <time.h>
//The following are used for fork().
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

//Load C++ libraries which are common to nearly all GAIA programs:
#include <iostream>
//The following libraries allow setw() to set widths in cout statements:
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
//Used for the random rand() function:
#include <cstdlib>
//Used for the min() function:
#include <algorithm>
//Used for string streams:
#include <sstream>

//Load custom library of LOTS of mathematical constants.
#include "double_constants.hpp"

//GSL libraries.
#ifdef GSL_HERE
  #include <gsl/gsl_linalg.h>
  #include <gsl/gsl_blas.h>
  #include <gsl/gsl_sf_erf.h>
  #include <gsl/gsl_sf_legendre.h>
  #include <gsl/gsl_errno.h>
  #include <gsl/gsl_fft_real.h>
  #include <gsl/gsl_fft_halfcomplex.h>
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_randist.h>
  #include <gsl/gsl_cdf.h>
#endif

//NetCDF libraries.
#ifdef NETCDF_HERE
  #include "netcdfcpp.h"
  static const int NC_ERR = 2;//Return this code to the OS in case of failure.
#endif

//For NFFT, NFFTLS.
#ifdef NFFT_HERE
  #include "io.h"
  #include "nfft.h"
  #include "utils.h"
  #include "ls.h"
  #include "usage.h"
#endif

struct input_options{
  string basefolder;//The folder (relative to GAIA/) that has input file.
  string filename;//input filename.
  string output_filename;//Occasionally, input data are saved using this.
  string output_folder;//Occasionally, input data are saved using this.
  int noload;//Set to 1 if the input isn't supposed to be loaded.
  int nocache;//Set to 1 if input.binindices isn't to be loaded from cache.
  int no_orbit_data;//Set to 1 if the input file has no orbit data.
  string cachefolder;//Folder where cached inputs are stored.
  long long segmentsize;//When loading input in segments, size of each segment.

  //Sometimes only some input data are used for inversion- other
  //points are reserved for verification.
  int split_input;
  long long chunk,buffer;//Used in split_input.

  //Constructor.
  input_options() : basefolder(),filename(),output_filename(),output_folder(),noload(),nocache(),no_orbit_data(),cachefolder(),segmentsize(),split_input(),chunk(),buffer() {}
};
typedef struct input_options input_options_s;

struct xy{
  //Each run of this program can produce multiple plots, so the data
  //are stored in a vector of vectors and the titles and units are also
  //stored in vectors.
  vector<string> titles;//Title on plots.
  vector<string> x_units,y_units;
  vector<int> log_x_axis,log_y_axis;//If 1, that plot uses log axis.
  vector< vector<string> > legends;//One for each plot's lines.
  //Vector of vectors of doubles for the output XY points.
  //1st index is the plot number- each plot is in a different file and pic.
  //2nd index is the dataseries- each plot can have multiple lines.
  //3rd index is the index of the x/y value in each dataseries.
  vector< vector< vector<double> > > x_values,y_values;
  //For labels, first applied to points in l-curve plots.
  //These have 1 fewer dimension than x_values because there's only
  //one set of labels per plot.
  vector< vector<double> > x_label_values,y_label_values;
  vector< vector<string> > labels;

  //Used to interpolate arbitrary xy timeseries to input time series.
  vector<long long> closest_indices;
  vector<double> differences;//for closest_indices, arbitrary - input.
  double numindices;
  double max_diff, min_diff;
  long long smallest_index, largest_index;
  vector<int> parameter_type;//copied from results.
  vector<double> parameter;//copied from results.
  long long current_parameter;
  string cachefolder;//Folder where preloaded arbitrary time series are stored.

  double xmin,xmax,ymin,ymax;

  //Constructor.
  xy() : titles(),x_units(),y_units(),log_x_axis(),log_y_axis(),legends(),x_values(),y_values(),x_label_values(),y_label_values(),labels(),closest_indices(),differences(),numindices(),max_diff(),min_diff(),smallest_index(),largest_index(),parameter_type(),parameter(),current_parameter(),cachefolder(),xmin(),xmax(),ymin(),ymax() {}
};
typedef struct xy xy_s;

struct topo{
  int type;//Controls how topography data are loaded.
  string basefolder;//Folder where topography data are stored.
  string filename;//Filename of topography data.
  vector<double> lat,lon,elev;

  //Used for SRTM30+ which is equal-grid and stored as 2-byte ints.
  vector< vector <short> > hr_elev;
  //hr_lats/lons are much shorter than hr_elev; they only store each row/column.
  vector<double> hr_lats, hr_lons;

  //1 - List of lats,lons,elevs.
  //2 - Elevs in 2D equal-grid grid.
  int format;

  //Constructor.
  topo() : type(),basefolder(),filename(),lat(),lon(),elev(),hr_elev(),hr_lats(),hr_lons(),format() {}
};
typedef struct topo topo_s;

struct love{
  int type;//Controls how Love numbers are loaded.
  string basefolder;//Folder where Love numbers are stored.
  string filename;//Filename of Love numbers.
  vector<double> h,k,l;//h,k,l Love numbers.

  //Constructor.
  love() : type(),basefolder(),filename(),h(),k(),l() {}
};
typedef struct love love_s;

struct temperatures{
  string basefolder;//Folder where temp data are stored.
  string filename;//Filename of temp data.
  vector<double> lat,lon,temp;

  //Constructor.
  temperatures() : basefolder(),filename(),lat(),lon(),temp() {}
};
typedef struct temperatures temp_s;

struct ancillary{
  topo_s topo;
  love_s love;
  temp_s temp;

  //Constructor.
  ancillary() : topo(),love(),temp() {}
};
typedef struct ancillary ancillary_s;

struct grid_options{
  int type;
  double lat_spacing, lon_spacing, cap_size;
  double boundaries[4];
  string cachefolder;//Folder where cached grids are stored.
  int noload;//Set to 1 if you don't want the grid loaded from cache.
  int area_type;//Specifies area of disks, squares, etc.

  int support_type;//Controls how supporting grid is formed.
  int support_lat_multiplier,support_lon_multiplier;

  //Needed because automatic cap sizing alters these.
  double original_lat_spacing, original_lon_spacing, original_cap_size;

  //Grid type 11, region type 3, analysis_list_choice 1 all use this point.
  double specified_lat,specified_lon;

  //In kilometers. When snapping grid/region to specified pt, the snap
  //is considered successful if the mismatch is less than this distance.
  double specified_tolerance;

  //Useful for cases when region code wraps around grid code.
  string base;//Usually either "grid" or "region".

  int topo_type;//Copy of ancillary.topo.type for grid cache substring.

  int latlon;//If 1, tells define_grid to make a latlon grid.
  double latlon_lat_multiplier, latlon_lon_multiplier;// >= 1, makes latlon denser than grid.
  int redundant_lats, redundant_lons;//0=avoid poles,lons0,360. 1=pts at poles AND lons 0,360 both.
  int global_latlon;//0=inherit unstructured grid boundaries, 1 = latlon grid is global.
  
  //Used only in region code.
  int parallel;//Switches on/off region parallel processing code.
  double overlap_distance;//Used to control how many "edge" points there are.

  int analysis_list_cache_only;//Skips any regions not in ls cache.
  int analysis_list_choice;//Used to determine which regions are analysed.
  //manual_list holds region indices that are used in some a_list_choices.
  vector<long long> manual_list;

  //Constructor.
  grid_options() : type(),lat_spacing(),lon_spacing(),cap_size(),boundaries(),cachefolder(),noload(),area_type(),support_type(),support_lat_multiplier(),support_lon_multiplier(),original_lat_spacing(),original_lon_spacing(),original_cap_size(),specified_lat(),specified_lon(),specified_tolerance(),base(),topo_type(),latlon(),latlon_lat_multiplier(),latlon_lon_multiplier(),redundant_lats(),redundant_lons(),global_latlon(),parallel(),overlap_distance(),analysis_list_cache_only(),analysis_list_choice(),manual_list() {}
};
typedef struct grid_options grid_options_s;

struct latlon{
  string filename;
  double mask;//Mask values occur on land.
  //Used for loading FES amplitudes and phases on 2D grid- obsolete/redundant?
  vector< vector<double> > amplitudes, phases;
  int current_index;//Used by read_raw_hybrid_cats, taken from model.constituent_indices.
  //Used for generic multi-parameter output on 2D grid.
  vector< vector< vector<double> > > outputs;
  vector<double> lat,lon;//Defines coords of 2D grid.
  vector< vector<double> > elev;//Elevation (m), probably that of nearby unstructured grid point.
  int discard_redundant_values;//Zero keeps them, anything else gets rid of redundant (lon=360) values in FES.

  //Records positions of original unstructured mascons.
  vector<double> mascon_lats,mascon_lons;

  //Constructor.
  latlon() : filename(),mask(),amplitudes(),phases(),current_index(),outputs(),lat(),lon(),elev(),discard_redundant_values(),mascon_lats(),mascon_lons() {}
};
typedef struct latlon latlon_s;

struct grid{
  //These variables are used by all grid types.
  vector<double> lat, lon, elev, x, y, z, distance_to_land, distance_to_ocean, wcover, wdepth;
  //Keeps a short, efficient list of the different lat bands in the grid.
  vector<double> different_lats;
  vector<long long> chosen_for_region;//debug!
  vector<double> area;//Area of each grid point for kg->cm conversion
  vector<double> temp;//Time-avg temp of each gr pt for kg->cm conversion(C)
  vector<double> density;//Density of each grid point for kg->cm conversion
  vector<double> cm2kg;//Used for cm-> kg conversion.
  vector<double> radius;//Radius(km) of each gr pt for disk mascons
  vector<double> radius_angle;//Radius(radians) for disk mascons

  //These variables are used to calculate the boundaries of bins
  //when the analysis method chosen uses bins (this is nearly always the case).
  double maxlat,minlat,maxlon,minlon;//Boundaries.
  double minlonsmall,maxlonsmall,minlonlarge,maxlonlarge;//Alt. boundaries.
  vector<double> lat_halfwidth;//Defined in terms of degrees or km.
  vector<double> lon_halfwidth;//Defined in terms of degrees or km.

  //Holds coordinates for supporting grid points. 1st index is main gr pt
  //that is "supported" by all the grid points in the 2nd index.
  //support_weights holds weighting factors for support grid points.
  vector<vector<double> > support_lat, support_lon, support_x, support_y, support_z, support_weights;

  //main_weights holds weighting factors for main grid points.
  vector<double> main_weights;

  //Every time the grid is rotated (in any direction) this variable
  //(in degrees) is used to perform the rotation. It's a temp var only.
  double rotate;

  //Holds copy of options so I don't have to pass grid_options around.
  grid_options_s options;
  
  //Holds equal-grid latlon grids.
  latlon_s latlon;

  //Only used by region code.

  //Contains all input data in this region.
  input_s input;

  //Each pt in regional_grid has an index to the original grid structure.
  //Some pts in regional_grid are "center" pts and some are "edge" pts.
  //First index is region #, second indices run over the grid points
  //in that region.
  vector< vector<long long> > all_grid_indices, center_grid_indices, edge_grid_indices;

  //Created by build_analysis_region_list().
  //Controls which regions are analyzed in main().
  vector<long long> indices_to_analyze;
  
  //Used to match FES pts with mascons.
  vector< vector<long long> > fes_lat_indices, fes_lon_indices;
  vector< vector< vector<long long> > > fes_index;
  vector<double> fes_areas;//Area for FES pts at different lats.

  //Constructor.
  grid() : lat(),lon(),elev(),x(),y(),z(),distance_to_land(),distance_to_ocean(),wcover(),wdepth(),different_lats(),chosen_for_region(),area(),temp(),density(),cm2kg(),radius(),radius_angle(),maxlat(),minlat(),maxlon(),minlon(),minlonsmall(),maxlonsmall(),minlonlarge(),maxlonlarge(),lat_halfwidth(),lon_halfwidth(),support_lat(),support_lon(),support_x(),support_y(),support_z(),support_weights(),main_weights(),rotate(),options(),latlon(),input(),all_grid_indices(),center_grid_indices(),edge_grid_indices(),indices_to_analyze(),fes_lat_indices(),fes_lon_indices(),fes_index(),fes_areas() {}
};
typedef struct grid grid_s;

struct period_info{
  int recognized;//Was current period recognized as a period defined here?
  int tide;//1 if this period is a tide (ex: M2) or not (ex: annual).
  int parameter_type;//Used in convert_title_to_parameter- for const, linear.
  string title;//Title of this period, if recognized.
  double period;//Period in seconds.
  double omega;//twoPi/(period in seconds) - plugs into cos(omega*t).
  int output_type;//Just like any entry in results.output_type.
  int sincos;//If period is encoded as sincos, this is 1. amp/phase = 0

  doodson_s doodson;//If tide, contains AT LEAST Doodson nums. Args/phase?

  //Constructor.
  period_info() : recognized(),tide(),parameter_type(),title(),period(),omega(),output_type(),sincos(),doodson() {}
};
typedef struct period_info period_info_s;

struct results{
  //The angular frequency ("omega") corresponding to every period is
  //stored here, in radians per second.
  vector<double> omegas;

  //Each run of this program can produce multiple plots, so the data
  //are stored in a vector of vectors and the titles and units are also
  //stored in vectors.
  vector<string> titles;//Title that ends up on plots.
  vector<string> units;//Units that are displayed on plots.
  string base_unit;//Usually "cm" but sometimes "nm/s^2"

  //Vector of vectors of doubles for mascon output.
  //1st index is the plot number- such as trend, period 1, period 2..
  //2nd index is the grid number.
  vector< vector<double> > outputs;

  //Vector of vectors of vectors of vectors of doubles for spherical harmonic output.
  //1st index is cycles through different parameters or monthly output fields.
  //2nd index is the degree "n", which starts at degree 0 so size() = n+1.
  //3rd index is the order "m".
  //4th index is 0 for Cnm, 1 for Snm.
  vector< vector< vector< vector<double> > > > harmonic_outputs;

  //debug - records how many times each entry in output is changed.
  vector<long long> output_changed;

  //convert_results_from_ampphase_to_sincos() sets this value based on the
  //original type of results.
  //1 - Results were already in sincos.
  //2 - Results were in amp/phase and were converted to sincos.
  int conversion_to_sincos;

  //convert_results_from_sincos_to_ampphase() sets this value based on the
  //original type of results.
  //1 - Results were already in ampphase.
  //2 - Results were in sincos and were converted to ampphase.
  int conversion_to_ampphase;

  //Holds copy of options so I don't have to pass analysis_options around.
  analysis_options_s options;

  //Holds output for 2D XY plots.
  xy_s xy;

  //Holds output for equal-grid latlon maps.
  latlon_s latlon;

  //Some maps (e.g. type 4) have labels for each pt.
  vector< vector<string> > labels;

  //Each map now has an optional error bar. Any errors < 0 are ignored.
  vector<double> error_bars;

  //Vector of vectors of doubles for markers on the plots.
  //1st index is the plot number- such as trend, period 1, period 2..
  //2nd index is the marker number in each plot.
  vector< vector<double> > marker_lats,marker_lons;

  //Used in minmax_output().
  double min,max;
  
  //Helpful for GMT and KMZ plots. Modified by is_polar().
  double minlat,maxlat,minlon,maxlon;

  vector<double> xyz;//For XYZ tristimulus values.
  vector<double> rgb;//For red,green,blue values. 0,0,0 = black; 1,1,1 = white.
  //Also for RGB maps:
  //0 - Obsolete. A single RGB scale, not uniformly normalized.
  //1 - Single RGB scale, uniformly normalized.
  //2 - Multiple RGB scales with different widths, all uniformly normalized.
  int rgb_choice;
  long long max_labels,max_widths;

  //Used when results_s holds a model:
  string title;//Title of this model, "FES 2004", "R. Ray", "Padman", etc.
  grid_s grid;//Allows model to be loaded onto its own unique grid.

  //1st index is M2, K1, etc. 2nd index is model.grid index.
  vector< vector<double> > amplitudes, phases, sin, cos;

  //The 1st index in amplitudes, phases refers to tide period. This tide
  //period is stored here, in seconds.
  vector<double> periods;

  string basefolder;//Folder (relative to GAIA/) that has model files.
  vector<string> filenames;//Model filenames.
  string cachefolder;//Folder where cached models are stored.

  //When comparing the model to tide gauge stations, this vector of vectors
  //holds the model's predicted sealevel height at each station.
  //1st index is the tide gauge station.
  //2nd index is the time index.
  vector< vector<double> > sealevels;//(Mean has been removed!)

  //Used if record_nearby_values==1.
  double avg_amp,avg_phase,closest_amp,closest_phase,closest_dist;

  vector<int> constituent_indices;//Records all indices that might be used by read_raw_hybrid_cats.

  //Constructor.
  results() : omegas(),titles(),units(),base_unit(),outputs(),harmonic_outputs(),output_changed(),conversion_to_sincos(),conversion_to_ampphase(),options(),xy(),latlon(),labels(),error_bars(),marker_lats(),marker_lons(),min(),max(),minlat(),maxlat(),minlon(),maxlon(),xyz(),rgb(),rgb_choice(2),max_labels(10),max_widths(3),title(),grid(),amplitudes(),phases(),sin(),cos(),periods(),basefolder(),filenames(),cachefolder(),sealevels(),avg_amp(),avg_phase(),closest_amp(),closest_phase(),closest_dist(),constituent_indices() {}
};
typedef struct results results_s;

struct plot_options{
  string outputfolder;//Folder where output is written.
  vector<string> output_files;//Filenames where outputs are written.
  string inputfolder;//Sometimes outputs are loaded, ex: load_output().

  int no_gmt_plots;//Set to 1 if you don't want GMT to make plots.
  string subfolder;//Folder where other output is written.
  
  //write_ascii_output creates a file with a list of filenames if this is the
  //final call of this fctn. This file lists the output filenames so
  //GMT knows which files to open.
  vector <string> just_the_filenames;

  string base;//Serves as the base for just_the_filenames.

  //Maps are made using this projection:
  int projection;

  int coastlines;

  //Plots with positive and negative values are scaled from +/- symmetric_limit.
  //Disable this feature by setting the limit <= 0.0.
  double symmetric_limit;

  //The number of digits in the color scale is controlled here:
  int scale_digits;

  //Determines color scheme for posters or for presentations.
  int color_scheme;

  //Left-justifies titles and adds (a), (b), etc., then runs montage.sh.
  int montage;

  //Sometimes the data range blurb isn't necessary.
  int blurb_disabled;

  //Phase maps will mask areas with amplitudes < error_bar in this way:
  int phase_mask;

  //Set to 1 to plot original mascon locations as brown points.
  int plot_mascons;

  int call_number;//Records which of cgaia.sh, c2gaia.sh, c3gaia.sh were used to run GAIA.
  
  //Constructor.
  plot_options() : outputfolder(),output_files(),inputfolder(),no_gmt_plots(),subfolder(),just_the_filenames(),base(),projection(),coastlines(),symmetric_limit(),scale_digits(),color_scheme(),montage(),blurb_disabled(),phase_mask(),plot_mascons(),call_number() {}
};
typedef struct plot_options plot_options_s;

//This is the include guard.

//Prototypes for the spectral color pipeline (Hughes and Williams 2010),
//excerpted from the full GAIA definitions.hpp.
void finish_flip_backgrounds(FILE *flip_fp);

void finish_trim(FILE *trim_fp);

string gaiafolder();

void match_x_axes(results_s &orig, results_s &match);

void interpolate(results_s &orig, results_s &match);

void closest_points(results_s &orig, results_s &match, long long i, long long j);

int is_this_point_closer(results_s &orig, results_s &match,
                         long long i, long long j, long long k);

void create_synthetic_plot(results_s &results, int choice, long long numpoints,
                           double x0, double delta, double param1, double param2);

void convert_timeseries_to_rgb(results_s &input, results_s &spectrum,
                               results_s &cie, plot_options_s &plot_options,
                               int choice, int nonmask_pts, int verbose);

void convert_timeseries_to_xyz(results_s &input, results_s &spectrum,
                               results_s &cie, plot_options_s &plot_options,
                               int choice, int nonmask_pts, int verbose);

void load_cie_functions(results_s &cie);

void spectrum2xyz(results_s &spectrum, results_s &cie, int verbose);

void xyz2rgb(results_s &spectrum);

void rescale_rgb(results_s &results);

void rescale_single_rgb(results_s &results, int verbose);

void rescale_single_xyz(results_s &results);

void gamma_correct(results_s &results);

void gamma_correct_single_rgb(results_s &results, int verbose);

void raise_y_to_power(results_s &results, double power);

void raise_single_y_to_power(results_s &results, double power);

void convert_spectrum_from_frequency_to_period(results_s &orig);

void spectra_init(results_s &spectra, grid_s &grid, int mod_choice, int &choice, int &norm_rgb);

void spectra_loop(results_s &timeseries, results_s &spectra, results_s &cie, results_s &coefs,
                  int &nonmask_pts, int mod_choice, int &choice, int &simple_detrend,
                  int &norm_rgb, double &parameter1, plot_options_s &plot_options);

void spectra_end(results_s &spectra, int mod_choice, int norm_rgb);

#endif
