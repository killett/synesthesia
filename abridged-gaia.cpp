#include "abridged-definitions.hpp"

void create_options(string &options_file, int argc, char *argv[]){
input_s input;
ancillary_s ancillary;
grid_options_s region_options;
grid_s grid;
analysis_options_s model_options;
results_s results;
plot_options_s plot_options;
insitu_s insitu;//For tide gauge data, etc.
input.options.filename="input_and_orbit_short3.txt";

input.options.split_input = 0;
input.options.buffer= 0;
input.options.chunk =  6681934;//Half of r5_0.05Hz_2002-2013_trimmed_50S

grid.options.type = 10;

grid.options.lat_spacing = 600;
grid.options.lon_spacing = grid.options.lat_spacing;

grid.options.cap_size = 1*deg2km;
grid.options.specified_lat = -77.36;//Location of 200km mascon which has anomalously high M2,S2 amps.
grid.options.specified_lon = -85.19;

//Global.
grid.options.boundaries[0] = -90.0;//min_lat
grid.options.boundaries[1] =  90.0;//max_lat
grid.options.boundaries[2] =   0.0;//min_lon
grid.options.boundaries[3] = 360.0;//max_lon */

double custom_lats[] = { 90.0 };
double custom_lons[] = { 180.0 };

grid.options.area_type = 2;

model_options.type = 0;

region_options.type = 2;

region_options.lat_spacing = 5000.0;
region_options.lon_spacing = region_options.lat_spacing;//Usually a good idea.

region_options.cap_size = 1700;

region_options.overlap_distance = 2050.0;

region_options.analysis_list_choice = 0;
region_options.specified_lat = grid.options.specified_lat;
region_options.specified_lon = grid.options.specified_lon;
long long manual_list[] = {1};//1st reg is 1, not 0.

results.options.type = 20001;

results.options.lambda = 20;

results.options.regparam1 = 1.0;
results.options.regparam2 = 48.0;

results.options.regparamint1 = 6;

results.options.time_resolution = 1;

results.options.convergence = 1E-15;

double parameter[] = {0,1,m2_period,k1_period,s2_period,o1_period,semiannual_period,annual_period};
long long p_type[] = {0,1,101,101,101,101,101,101};
double lambdas[] = {20,20,18.5,65,550,47.5,20,20};// */

results.options.h_matrix_choice = 202;

grid.options.support_type=0;

results.options.zone_of_influence_choice = 1;
results.options.zone_radius = region_options.overlap_distance;

insitu.comparison_type = 0;

results.options.upward_continuation_verification_choice = 0;
results.options.skip_upward_continuation_during_verification = 0;
results.options.folder_for_skipping = "20111016-165723/";

grid.options.redundant_lats=1;
grid.options.redundant_lons=1;

grid.options.global_latlon=0;

grid.options.latlon_lat_multiplier=1.0;//Crude but fast, nmax<=200. 
grid.options.latlon_lon_multiplier=grid.options.latlon_lat_multiplier;

plot_options.no_gmt_plots = 0;

//Global maps will use this projection:
// 1 - Robinson - A projection designed to be aesthetically appealing.
// 2 - Winkel Tripel - Minimizes sum of area, angle, and distance distortions.
// 3 - Mollweide - Pseudocylindrical equal area projection.
// 4 - Miller - Cylindrical projection, similar to Mercator but extends to the poles. 
//(Polar maps automatically use the azimuthal equidistant projection.)
plot_options.projection = 2;

//Coastlines:
// 1 - GMT pscoast $pscoast_res
// 2 - GMT pscoast $pscoast_res + Rignot2011 InSAR.
//     Rignot, E., J. Mouginot, and B. Scheuchl. 2011. Antarctic Grounding Line Mapping from Differential Satellite Radar Interferometry, Geophyical Research Letters, 38, L10504, doi:10.1029/2011GL047109.
//     Rignot, E., J. Mouginot, and B. Scheuchl. 2011. MEaSUREs Antarctic Grounding Line from Differential Satellite Radar Interferometry. Boulder, Colorado USA: NASA EOSDIS DAAC at NSIDC. http://nsidc.org/data/nsidc-0498.htm
plot_options.coastlines = 1;

//Plots with positive and negative values are scaled from +/- symmetric_limit.
//Disable this feature by setting the limit <= 0.0.
plot_options.symmetric_limit = -1.0;

//The color scale will display this many digits after the decimal place.
//If the color scale max/min are too big/small to display
//with the specified digits, scientific notation will be used.
//If scale_digits==1(or "n"), scientific notation will be used if:
//the color scale max or min are >= 1000.0 (or 10^(n+2))
//the color scale max or min are < 0.1  (or 10^(-n))
//Otherwise, fixed notation will be used.
//If abs value of max/min of data less than num of digits, GMT fails to plot!
plot_options.scale_digits = 4;

//All plots will use this color scheme:
// 1 - POSTER Black text, page is white.
// 2 - White text, page is black.
plot_options.color_scheme = 1;

//Setting this to 1 left-justifies titles, adds (a)'s, and runs montage.sh:
plot_options.montage = 0;

//Phase maps will mask areas with amplitudes < error_bar in this way:
// 1 - POSTER Low amplitude areas are abruptly masked with off-white.
// 2 - Low amplitude areas are abruptly masked with dark gray.
//11 - POSTER Low amplitude areas gradually fade to white.
//12 - Low amplitude areas are gradually dimmed, approaching black.
plot_options.phase_mask = 12;

//Set to 1 to plot original mascon locations as brown points.
plot_options.plot_mascons = 1;

int main(int argc, char *argv[]){
  long long i,k;//Counters.
  string options_file(gaiafolder());
  input_s input, unused_input;
  ancillary_s ancillary;
  grid_s grid, region;
  results_s model, results, regional_results,overlap_results;
  plot_options_s plot_options,temp_options,overlap_plot_options;
  insitu_s insitu;//For tide gauge data, etc.
  //These variables are used to display run time of program at end.
  time_t start_time,end_time;
  timespan_conversion_s timespan_conversion;//From secs to human units

  char s[max_length];

  string temp_string;

  //Used to search the overlap folder for filenames of each region's files.
  DIR *cacheDIR;
  struct dirent *cacheDirEnt;
  string newfile;

  start_time = time(NULL);//Record time at start of program.

  create_options(options_file, argc, argv);

  //Load options from file produced by create_options.txt
  load_options(options_file, input, region.options, ancillary, grid, model.options, results, plot_options, insitu,1);
  cout<<"call number after load_options: "<<results.options.call_number<<endl;

  cout<<"Saving output to "<<plot_options.outputfolder<<endl;

  load_ancillary(ancillary, 1);

  define_grid(grid, ancillary, results);

  //If requested, build latlon grid of "sufficient" density, copy output to it.
  copy_results_to_latlon(results, grid);

  //Save latlon output to disk as ascii files.
  plot_options.base = "latlon_data";
  write_structured_output(results,plot_options);
  //Save latlon output to disk as NetCDF files, and write GMT,KMZ scripts.
  if(write_netcdf_output(results,grid,plot_options,1) != 0) cout<<"!!!WARNING!!! Problem writing NetCDF output files!"<<endl;

  end_time = time(NULL);//Record time at end of program.
  timespan_conversion.seconds = (long long)(end_time - start_time);
  cout<<"GAIA v"<<fixed<<setprecision(3)<<GAIA_VERSION<<" finished after"<<sec2human(timespan_conversion)<<"."<<endl;

  return 0;
}
