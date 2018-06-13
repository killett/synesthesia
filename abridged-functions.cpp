#include "abridged-definitions.hpp"

void write_gmt_cpt(results_s &results,
                   plot_options_s &plot_options,
                   FILE *new_fp,
                   long long i){
  /**********************************************************************
  Purpose: This function writes the GMT command for creating the cpt file
          used for making the color scale.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          new_fp - file pointer to current gmt script file
          i - (output) index of parameter being scaled
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  int numlevels = 10;
  string palette_guide("#0-360phase=cyclic,amps=wysiwyg,trend=polar,topo=relief,pts=rainbow");
  string overflow_guide("#-E=both,-Ef=top,-Eb=bottom");

  //Search for min/max, store in results.min/max.
  minmax_output(results,i,0);//Map index, quiet.
  
  if(results.options.output_choice == 1){
    if(results.options.output_type.at(i)==104){
      fprintf(new_fp,"palette=\" -Ccyclic \" %s\n",palette_guide.c_str());
      fprintf(new_fp,"upper_limit=\"%.6f\"\n",180.0);
      fprintf(new_fp,"lower_limit=\"%.6f\"\n",-180.0);
      //Don't need decimal places to display these integer limits.
      //BUT IT ALSO CHANGES EVERY MAP AFTER THIS. FIX THAT!
      //fprintf(new_fp,"scale_format=\" --D_FORMAT=%%.0f \"\n");
    }
    else{
      fprintf(new_fp,"palette=\" -Crainbow \" %s\n",palette_guide.c_str());
      fprintf(new_fp,"upper_limit=\"%.6f\"\n",results.max);
      fprintf(new_fp,"lower_limit=\"%.6f\"\n",results.min);
    }
    fprintf(new_fp,"overflow=\"\" %s\n",overflow_guide.c_str());
    fprintf(new_fp,"numlevels=\"%d\"\n",numlevels);
    fprintf(new_fp,"# The \"| bc -l\" helps bash deal with floating point numbers.\n");
    fprintf(new_fp,"interval=\"$(echo \"($upper_limit - $lower_limit) / $numlevels\" | bc -l)\"\n");
    fprintf(new_fp,"$gmt_prefix makecpt $palette -T$lower_limit/$upper_limit/$interval -Z > map.cpt\n");
  }
  else if(results.options.output_choice == 5){
    //If this is a phase, make a cyclic, amplitude-shaded cpt file.
    if(results.options.output_type.at(i)==104){
      cout<<"NOTE: Phases should only use the cyclic colorscale if the min/max span [-180,180] or [0,360]!"<<endl;
      fprintf(new_fp,"palette=\" -Ccyclic -I \" %s\n",palette_guide.c_str());
      fprintf(new_fp,"limits=\"\" #Default limits include all values (not just in subset), aren't always symmetric.\n");
      fprintf(new_fp,"overflow=\"\" %s\n",overflow_guide.c_str());
      numlevels = 13;//If 360 degrees, is divided by 12 (+1 for value 0.0 because 180 is even) nicely
      fprintf(new_fp,"numlevels=\"%d\"\n",numlevels);
      fprintf(new_fp,"$gmt_prefix grd2cpt $data_name $palette $limits -E$numlevels -Z > map.cpt\n");
    }
    else if(results.rgb.size() != 3){
      fprintf(new_fp,"subset_name=$data_name\"_subset\" #Weird to append, but works with ../pl..\n");
      fprintf(new_fp,"$gmt_prefix grdcut $data_name -G$subset_name $actual_range #Only use subset values for scale.\n");
      fprintf(new_fp,"$gmt_prefix grdinfo -L0 $subset_name > subset_grdinfo #Extract data range from subset.\n");
      fprintf(new_fp,"data_min_e=$(awk '/z_min: /{printf \"%%.'$digits'e\\n\", $3}' subset_grdinfo)\n");
      fprintf(new_fp,"data_max_e=$(awk '/z_max: /{printf \"%%.'$digits'e\\n\", $5}' subset_grdinfo)\n");
      fprintf(new_fp,"data_min_f=$(awk '/z_min: /{printf \"%%.'$digits'f\\n\", $3}' subset_grdinfo)\n");
      fprintf(new_fp,"data_max_f=$(awk '/z_max: /{printf \"%%.'$digits'f\\n\", $5}' subset_grdinfo)\n");
      fprintf(new_fp,"notation=\"f\" #By default, numbers appear in floating point format.\n");
      fprintf(new_fp,"data_min_print=$data_min_f #By default, numbers appear in floating point format.\n");
      fprintf(new_fp,"data_max_print=$data_max_f #By default, numbers appear in floating point format.\n");
      fprintf(new_fp,"palette=\" -Cwysiwyg \" %s\n",palette_guide.c_str());
      fprintf(new_fp,"symmetric_limit=\"%.6f\"\n",plot_options.symmetric_limit);
      fprintf(new_fp,"if [[ $(bc <<< \"$symmetric_limit <= 0.0\") == 1 || $(bc <<< \"$data_min_f >= 0.0\") == 1 || $(bc <<< \"$data_max_f <= 0.0 \") == 1 ]]\n");
      fprintf(new_fp,"then\n");
      fprintf(new_fp,"  upper_limit=$data_max_f\n");
      fprintf(new_fp,"  lower_limit=$data_min_f\n");
      fprintf(new_fp,"else\n");
      fprintf(new_fp,"  upper_limit=$symmetric_limit\n");
      fprintf(new_fp,"  lower_limit=-$upper_limit\n");
      fprintf(new_fp,"fi\n");
      fprintf(new_fp,"#upper_limit=0.0\n");
      fprintf(new_fp,"#lower_limit=0.0\n");
      fprintf(new_fp,"limits=\" -L$lower_limit/$upper_limit \"\n");
      fprintf(new_fp,"#limits=\"\" #Default limits include all values (not just in subset), aren't always symmetric.\n");
      fprintf(new_fp,"#limits=\" -T= \" #Symmetric limits that include all values (not just in subset).\n");
      fprintf(new_fp,". ./overflow.sh\n");
      fprintf(new_fp,"if [[ $(bc <<< \"$data_min_f >= 0.0\") == 1 ]]\n");
      fprintf(new_fp,"then\n");
      fprintf(new_fp,"  #If all values are non-negative, they might cluster near bottom of\n");
      fprintf(new_fp,"  #scale which is dark so use off-white coastlines.\n");
      fprintf(new_fp,"  coast_color=\"gray82\"\n");
      fprintf(new_fp,"else\n");
      fprintf(new_fp,"  #Otherwise, they might cluster near the scale center which is light so use dark coastlines.\n");
      fprintf(new_fp,"  coast_color=\"gray10\"\n");
      fprintf(new_fp,"fi\n");
      fprintf(new_fp,"numlevels=\"%d\"\n",numlevels);
      fprintf(new_fp,"$gmt_prefix grd2cpt $data_name $palette $limits -E$numlevels -Z > map.cpt\n");
    }
    else fprintf(new_fp,"  coast_color=\"gray82\"\n");//For RGB maps.
  }
  else cout<<"!!!!WARNING!!!!!! results.options.output_choice "<<results.options.output_choice<<" isn't recognized."<<endl;
  plot_options.blurb_disabled=0;//Changed my mind- better to always see it.
}

void write_gmt_defs(FILE *new_fp){
  /**********************************************************************
  Purpose: This function writes some clarifying definitions that help
          me to consistently write working GMT data plotting commands.
          2017-04-04 Update: Also added GMT4/5 compatibility which uses
                             the "GMT5" define from definitions.hpp.
                             THIS ISN'T FINISHED YET!
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          new_fp - file pointer to current gmt script file
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  fprintf(new_fp,"#######################################################\n");
  fprintf(new_fp,"#Clarifying definitions. Do not change!################\n");
  fprintf(new_fp,"start=\" -K \" #Should always redirect using > to write new PS.\n");
  fprintf(new_fp,"middle=\" -O -K \" #Should always redirect using >> to append to PS.\n");
  fprintf(new_fp,"end=\" -O \" #Should always redirect using >> to append to PS.\n");
  fprintf(new_fp,"#######################################################\n");
  //GMT5 requires all GMT4 commands to be preceded by "gmt " which this
  //simple change addresses, but it doesn't address all the other changes in GMT5's syntax.
  #ifdef GMT5
    fprintf(new_fp,"gmt5=1 #1/0 = GMT v5/4. GMTv5 support not finished.\n");
  #else
    fprintf(new_fp,"gmt5=0 #1/0 = GMT v5/4. GMTv5 support not finished.\n");
  #endif
  fprintf(new_fp,"if [ $gmt5 != 0 ]\n");
  fprintf(new_fp,"then\n");
  fprintf(new_fp,"  gmt_prefix=\"gmt \"\n");
  fprintf(new_fp,"else\n");
  fprintf(new_fp,"  gmt_prefix=\"\"\n");
  fprintf(new_fp,"fi\n");
}

void write_gmt_map_data(results_s &results,
                        grid_s &grid,
                        plot_options_s &plot_options,
                        FILE *new_fp,
                        string title,
                        long long i){
  /**********************************************************************
  Purpose: This function writes the GMT data plotting command.
          If the title is blank, this is assumed to be a KMZ map so
          the map fills the page in the cylindrical equidistant
          projection that works in Google Earth.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          new_fp - file pointer to current gmt script file
          title - string to be plotting at top of map data. If blank,
                  title is suppressed and Google Earth KMZ is output.
          i - (output) index of parameter being mapped.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  //Is this a polar or global plot?
  int polar = is_polar(results,grid);//0 - global, 1 - north pole, 2 - south pole.

  if(i==0){
    fprintf(new_fp,"#Set resolution, coast_file, coast_thickness, and coastlines\n");
    fprintf(new_fp,"#on first map only because they should be universal.\n");
    string coast_file(gaiafolder());
    coast_file.append("data/ancillary/Rignot/InSAR_GL_Antarctica.txt");
    fprintf(new_fp,"coast_file=\"%s\"\n",coast_file.c_str());
    //Only latlon data determines resolution using delta_lat.
    if(results.options.output_choice == 5){
      double delta_lat = fabs(results.latlon.lat.at(1) - results.latlon.lat.at(0));
      if(delta_lat < 0.4){//Latlon lat spacing controls the resolution.
        fprintf(new_fp,"resolution=\" -E50 \" #50/2000 is low/high quality.\n");
        fprintf(new_fp,"pscoast_res=\" -Df+ \"\n");
      }
      else if(delta_lat < 0.9){//Latlon lat spacing controls the resolution.
        fprintf(new_fp,"resolution=\" -E50 \" #50/2000 is low/high quality.\n");
        fprintf(new_fp,"pscoast_res=\" -Df+ \"\n");
      }
      else{
        fprintf(new_fp,"resolution=\" -E50 \" #50/2000 is low/high quality.\n");
        fprintf(new_fp,"pscoast_res=\" -Di+ \"\n");
      }
    }
    else if(results.options.output_choice == 1 or results.options.output_choice == 4){
      fprintf(new_fp,"resolution=\" -E50 \" #50/2000 is low/high quality.\n");
      fprintf(new_fp,"pscoast_res=\" -Di+ \"\n");
    }
    else cout<<"!!!!WARNING!!!!!! results.options.output_choice "<<results.options.output_choice<<" isn't recognized."<<endl;
    fprintf(new_fp,"pscoast_res_orig=$pscoast_res #Don't want USA maps to repeatedly add -N2.\n");
    fprintf(new_fp,"pscoast_thk=\"0.6\"\n");
    fprintf(new_fp,"coast_thk=\"0.009\"\n");
    if(polar==1) plot_options.coastlines = 1;//InSAR is only in Antarctica, so disable for NP plots.
    fprintf(new_fp,"coastlines=%d #1:pscoast, 2:pscoast+InSAR.\n",plot_options.coastlines);
  }
  fprintf(new_fp,"#coast_color is gray82 for off-white, or gray10 for dark coastlines.\n");
  
  //Adjust max/min latitudes for mapping points, otherwise points on edge aren't visible.
  if(results.options.output_choice == 1 or results.options.output_choice == 4){
    double buffer = 5;
    if(results.maxlat <= 90-buffer) results.maxlat += buffer;
    if(results.minlat >= -90+buffer) results.minlat -= buffer;
  }

  //Record text formats, which are the same for all projections and data types.
  fprintf(new_fp,"title_format=\"0 0 30 0 0 MC\"\n");
  fprintf(new_fp,"blurb_format=\"0 0 15 0 1 ML\"\n");
  fprintf(new_fp,"units_format=\"0 0 13 0 0 MC\"\n");
  //Record units for the scale, which are the same for all projections and data types.
  fprintf(new_fp,"scale_units=\"%s\"\n",results.units.at(i).c_str());

  fprintf(new_fp,"misc_range=\" -R0/1/0/1 -JX1c \"\n");
  fprintf(new_fp,"#grdcut requires actual limits, but if grdimage uses them: GMT Fatal Error: grdimage could not allocate memory [21.69 Gb, n_items = 5823567396]\n");
  fprintf(new_fp,"minlon=%.3f\n",0.0);//results.minlon);
  fprintf(new_fp,"maxlon=%.3f\n",360.0);//results.maxlon);
  fprintf(new_fp,"minlat=%.3f\n",-90.0);//results.minlat);
  fprintf(new_fp,"maxlat=%.3f\n",90.0);//results.maxlat);

  if(!title.empty()){
    //All projections are always available.
    if(i==0){//Only print this guide for the first map.
      fprintf(new_fp,"#Global projections:\n");
      fprintf(new_fp,"#    1 - Robinson\n");
      fprintf(new_fp,"#    2 - Winkel Tripel\n");
      fprintf(new_fp,"#    3 - Mollweide\n");
      fprintf(new_fp,"#    4 - Miller\n");
      fprintf(new_fp,"#Polar projections:\n");
      fprintf(new_fp,"#  101 - N. Azimuthal Equidistant\n");
      fprintf(new_fp,"#  102 - S. Azimuthal Equidistant\n");
      fprintf(new_fp,"#Specific regions:\n");
      fprintf(new_fp,"# 1001 - North America\n");
      fprintf(new_fp,"# 1002 - South America\n");
      fprintf(new_fp,"# 1003 - Africa\n");
      fprintf(new_fp,"# 1004 - Greenland\n");
      fprintf(new_fp,"# 1005 - South Asia\n");
      fprintf(new_fp,"# 1006 - Australia\n");
      fprintf(new_fp,"# 1007 - Europe\n");
      fprintf(new_fp,"# 1101 - Contiguous United States\n");
      fprintf(new_fp,"# 1102 - California\n");
    }
    if(polar==0){
      if(i>0) fprintf(new_fp,"#");
      fprintf(new_fp,"projection_choice=%d\n",plot_options.projection);
    }
    else{
      if(polar==1){//North pole.
        if(i>0) fprintf(new_fp,"#");
        fprintf(new_fp,"projection_choice=101\n");
      }
      else{//South pole.
        if(i>0) fprintf(new_fp,"#");
        fprintf(new_fp,"projection_choice=102\n");
      }
    }
    fprintf(new_fp,"standard_circle=0 #1=all specific regions use standard circular projection.\n");
    fprintf(new_fp,"standard_rect=0 #1=all specific regions use standard rectangular projection.\n");
    fprintf(new_fp,". ./projections.sh\n");
    fprintf(new_fp,"if [ $projection_choice == 101 ]\n");
    fprintf(new_fp,"then\n");
    if(polar==1) fprintf(new_fp,"  minlat=%.3f\n",results.minlat);
    else         fprintf(new_fp,"  minlat=%.3f\n",0.0);
    fprintf(new_fp,"  actual_range=\" -R0.0/360.0/$minlat/90.0 \"\n");
    fprintf(new_fp,"  polar_radius=$(bc <<< \"scale=5; 90-$minlat\")\n");
    fprintf(new_fp,"  projection=\" -JE0/90.0/$polar_radius/${map_width}c \" #N. Azimuthal Equidistant\n");
    fprintf(new_fp,"elif [ $projection_choice == 102 ]\n");
    fprintf(new_fp,"then\n");
    if(polar==2) fprintf(new_fp,"  maxlat=%.3f\n",results.maxlat);
    else         fprintf(new_fp,"  maxlat=%.3f\n",0.0);
    fprintf(new_fp,"  actual_range=\" -R0.0/360.0/-90.0/$maxlat \"\n");
    fprintf(new_fp,"  polar_radius=$(bc <<< \"scale=5; 90+$maxlat\")\n");
    fprintf(new_fp,"  projection=\" -JE0/-90.0/$polar_radius/${map_width}c \" #S. Azimuthal Equidistant\n");
    fprintf(new_fp,"fi\n");
    fprintf(new_fp,"range=\" -R$minlon/$maxlon/$minlat/$maxlat \"\n");
    fprintf(new_fp,"map_pos=\" -Xa${map_x}c -Ya${map_y}c \"\n");
    if(results.rgb.size() == 3 and results.rgb_choice >= 2) fprintf(new_fp,"scale_width=1.2 #Override for RGB maps.\n");
    fprintf(new_fp,"scale_pos=\" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c \"\n");
    fprintf(new_fp,"units_x=$(bc <<< \"scale=5; $scale_x+$scale_width/2\")\n");
    fprintf(new_fp,"units_y=$(bc <<< \"scale=5; $scale_y+$scale_length/2+0.8\")\n");
    fprintf(new_fp,"units_pos=\" -Xa${units_x}c -Ya${units_y}c \"\n");
    fprintf(new_fp,"blurb_pos=\" -Xa${blurb_x}c -Ya${blurbs_y}c \"\n");
    fprintf(new_fp,"blurb2_pos=\" -Xa${blurb2_x}c -Ya${blurbs_y}c \"\n");
  }
  else{
    fprintf(new_fp,"range=\" -R0/360/%.3f/%.3f \"\n",results.minlat,results.maxlat);
    fprintf(new_fp,"actual_range=\" -R%.3f/%.3f/%.3f/%.3f \"\n",results.minlon,results.maxlon,results.minlat,results.maxlat);
    fprintf(new_fp,"projection=\" -JQ0/29c \" #Cylindrical Equidistant\n");
    fprintf(new_fp,"map_pos=\" -Xa0.1c -Ya5c \"\n");
    fprintf(new_fp,"scale_y=10.0\n");
    fprintf(new_fp,"scale_length=10.0\n");
    if(results.rgb.size() != 3 or results.rgb_choice < 2) fprintf(new_fp,"scale_width=0.6\n");
    else fprintf(new_fp,"scale_width=1.2\n");
    fprintf(new_fp,"scale_x=5.0\n");
    fprintf(new_fp,"scale_pos=\" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c \"\n");
    fprintf(new_fp,"units_x=$(bc <<< \"scale=5; $scale_x+$scale_width/2\")\n");
    fprintf(new_fp,"units_y=$(bc <<< \"scale=5; $scale_y+$scale_length/2+0.8\")\n");
    fprintf(new_fp,"units_pos=\" -Xa${units_x}c -Ya${units_y}c \"\n");
    fprintf(new_fp,"title_pos=\" -Xa-4c -Ya-5c \"\n");
    fprintf(new_fp,"blurb_pos=\" -Xa1c -Ya5c \"\n");
  }

  //Create cpt file for coloring map.
  write_gmt_cpt(results,plot_options,new_fp,i);

  if(results.options.output_choice == 5){
    //If this is a phase, make a new grd file that can shade
    //the phases using the value of the amplitudes, which MUST
    //be the file immediately preceding this one.
    if(results.options.output_type.at(i)==104){
      //Search for min/max of previous amplitude map, store in results.min/max.
      minmax_output(results,i-1,0);//Map index, quiet.
      //Should be results.error_bar(i-1);
      double masked_amp = results.min + 0.1*(results.max-results.min);
      fprintf(new_fp,"masked_amp=%.2f\n",masked_amp);
      plot_options.phase_mask = 12;
      cout<<"NOTE: phase_mask was forced to "<<plot_options.phase_mask<<endl;
      if(plot_options.phase_mask == 1){//POSTER Low amplitude areas are abruptly masked with off-white.
      }
      else if(plot_options.phase_mask == 2){//Low amplitude areas are abruptly masked with dark gray.
      }
      if(plot_options.phase_mask == 11){//POSTER Low amplitude areas gradually fade to white.
        //ATAN function maps positive values to [0,pi/2] so reverse and map to [0,1]
        fprintf(new_fp,"$gmt_prefix grdmath $previous_data_name $masked_amp GT ATAN 1.5708 DIV -2 MUL 1 ADD = amp_shading.grd\n");//i-1 is the amplitude data file.
        //Plot the shading file just to see what it looks like.
        fprintf(new_fp,"$gmt_prefix grd2cpt amp_shading.grd -Cgray -Z > shading.cpt\n");
        //KMZ scripts don't define $boundary, so it's not used here.
        fprintf(new_fp,"$gmt_prefix grdimage amp_shading.grd -B:.\"Amplitude shading\": $resolution $range $projection $map_pos -Cshading.cpt $start > amp_shading.ps\n");
        fprintf(new_fp,"$gmt_prefix psscale -Cshading.cpt -B/:#: -L --D_FORMAT=%%.2e $overflow $scale_pos -A $end >> amp_shading.ps\n");
        fprintf(new_fp,"$gmt_prefix grd2cpt $data_name $palette -E$numlevels -I -Z > map.cpt\n");
        fprintf(new_fp,"map_pos=$map_pos\" -Iamp_shading.grd \"\n");
        //Record lower amplitude used for masking phase values.
        fprintf(new_fp,"blurb_contents=\"Faded for amp. < %.2f %s\"\n",masked_amp,results.units.at(i-1).c_str());
        fprintf(new_fp,"coast_color=\"gray10\"\n");//Dark gray coastlines work well for phases that fade to white.
      }
      else if(plot_options.phase_mask == 12){//Low amplitude areas are gradually dimmed, approaching black.
        //ATAN function maps positive values to [0,pi/2] so map to [-1,0]
        fprintf(new_fp,"$gmt_prefix grdmath $previous_data_name $masked_amp GT ATAN 1.5708 DIV 2 MUL 1 SUB = amp_shading.grd\n");//i-1 is the amplitude data file.
        //Plot the shading file just to see what it looks like.
        fprintf(new_fp,"$gmt_prefix grd2cpt amp_shading.grd -Cgray -Z > shading.cpt\n");
        //KMZ scripts don't define $boundary, so it's not used here.
        fprintf(new_fp,"$gmt_prefix grdimage amp_shading.grd -B:.\"Amplitude shading\": $resolution $range $projection $map_pos -Cshading.cpt $start > amp_shading.ps\n");
        fprintf(new_fp,"$gmt_prefix psscale -Cshading.cpt -B/:#: -L --D_FORMAT=%%.2e $overflow $scale_pos -A $end >> amp_shading.ps\n");
        fprintf(new_fp,"$gmt_prefix grd2cpt $data_name $palette -E$numlevels -I -Z > map.cpt\n");
        fprintf(new_fp,"map_pos=$map_pos\" -Iamp_shading.grd \"\n");
        //Record lower amplitude used for masking phase values.
        fprintf(new_fp,"blurb_contents=\"Dimmed for amp. < $masked_amp %s\"\n",results.units.at(i-1).c_str());
        fprintf(new_fp,"coast_color=\"gray82\"\n");//Off-white coastlines work well for phases that are dimmed, approaching black.
      }
      else cout<<"!!!WARNING!!! The phase_mask "<<plot_options.phase_mask<<" wasn't recognized!"<<endl;
      fprintf(new_fp,"scale_format=\" --D_FORMAT=%%.0f \" #Full 0-360 scale is best w/ ints.\n");
    }
    else{
      if(results.rgb.size() == 3) fprintf(new_fp,"#");
      fprintf(new_fp,". ./notation.sh\n");
      fprintf(new_fp,"blurb_contents=\"Data range: [$data_min_print, $data_max_print] $scale_units\"\n");
    }
  }

  if(!title.empty()){
    if(results.options.output_choice == 1){
      fprintf(new_fp,"title=\"%s\"\n",title.c_str());
      //Create basemap, with title on top.
      fprintf(new_fp,"$gmt_prefix psbasemap $boundary:.\"$title\": $range $projection $map_pos $start > $plot_base.ps\n");
      write_gmt_coastlines(new_fp);
      fprintf(new_fp,"sym_size=\"0.3c\"\n");
      fprintf(new_fp,"$gmt_prefix psxy -N $data_name -bclongitudes/latitudes/output -Sc$sym_size -Cmap.cpt $range $projection $map_pos $middle >> $plot_base.ps\n");
      fprintf(new_fp,"if [ $montage != 0 ]\n");
      fprintf(new_fp,"then\n");
      fprintf(new_fp,"  title=${prefixes[$index]}\" \"$title\n");
      fprintf(new_fp,"  title_format=\"0 0 30 0 0 ML\" #Left-justify so montage titles are uniform.\n");
      fprintf(new_fp,"  title_x=$(bc <<< \"scale=5; $blurb_x-0.1\")\n");
      fprintf(new_fp,"else\n");
      fprintf(new_fp,"  title_x=$(bc <<< \"scale=5; $map_x+$map_width/2\")\n");
      fprintf(new_fp,"fi\n");
      fprintf(new_fp,"title_pos=\" -Xa${title_x}c -Ya${title_y}c \"\n");
      fprintf(new_fp,"echo $title_format $title | $gmt_prefix pstext -N $title_pos $misc_range $middle >> $plot_base.ps\n");
      
      //Draw color scale with units printed above.
      write_gmt_colorscale(new_fp,results,0);//0 = next to map for GMT PDF output.

      //Print blurb about data range, or masked amplitudes for phase plots.
      if(plot_options.blurb_disabled) fprintf(new_fp,"blurb_contents=\"\"\n");
      fprintf(new_fp,"echo $blurb_format $blurb_contents | $gmt_prefix pstext -N $blurb_pos $misc_range $middle >> $plot_base.ps\n");
      
      //Only print error bars if they're the right length, and this one is > 0.0.
      //Need to have separate copies for unstructured and latlon maps because
      //in each case the "right length" is defined by a different multivector.
      int blurb2_written=0;
      if(results.outputs.size() == results.error_bars.size() and !results.outputs.empty()){
        if(results.error_bars.at(i) > 0.0){
          blurb2_written=1;
          fprintf(new_fp,"blurb2_contents=\"Error bar: %.1f $scale_units\"\n",results.error_bars.at(i));
        }
      }
      if(!blurb2_written) fprintf(new_fp,"blurb2_contents=\"\" #Error bar: N/A $scale_units\"\n");
      fprintf(new_fp,"echo $blurb_format $blurb2_contents | $gmt_prefix pstext -N $blurb2_pos $misc_range $end >> $plot_base.ps\n");
    }
    else if(results.options.output_choice == 4){
      cout<<"!!!WARNING!!! 2D PLOTS SHOULD GO TO WRITE_GMT_PLOT_DATA!"<<endl;
    }
    else if(results.options.output_choice == 5){
      //Plot data, with title on top.
      fprintf(new_fp,"title=\"%s\"\n",title.c_str());
      if(results.rgb.size()==3 and results.latlon.outputs.size()==3) fprintf(new_fp,"$gmt_prefix grdimage red.nc green.nc blue.nc $boundary $resolution $range $projection $map_pos $start > $plot_base.ps\n");
      else fprintf(new_fp,"$gmt_prefix grdimage $data_name $boundary $resolution $range $projection $map_pos -Cmap.cpt $start > $plot_base.ps\n");
      write_gmt_coastlines(new_fp);
      //If marker_lats multivectors have the right size, plot markers in the color "sandy brown".
      //Plot markers before mascons because mascons are smaller than markers.
      if(results.marker_lats.size() == results.latlon.outputs.size() and results.marker_lons.size() == results.latlon.outputs.size() and !results.latlon.outputs.empty()){
        if(!results.marker_lats.at(i).empty() and !results.marker_lons.at(i).empty()){
          fprintf(new_fp,"$gmt_prefix psxy -N $data_name -bcmarker_lons/marker_lats -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n");
        }
      }
      fprintf(new_fp,"#Uncomment to put a marker at echoed coords, given as lon lat:\n");
      fprintf(new_fp,"#echo -85.19 -77.36 | $gmt_prefix psxy -N -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n");
      //If requested and mascon_lats/lons vectors aren't empty, plot mascon centers in the color "saddle brown".
      if(plot_options.plot_mascons != 0 and !results.latlon.mascon_lats.empty()) fprintf(new_fp,"$gmt_prefix psxy $data_name -bcmascon_lons/mascon_lats -Sc0.01c -G139/69/19 $range $projection $map_pos $middle >> $plot_base.ps\n");
      fprintf(new_fp,"if [ $montage != 0 ]\n");
      fprintf(new_fp,"then\n");
      fprintf(new_fp,"  title=${prefixes[$index]}\" \"$title\n");
      fprintf(new_fp,"  title_format=\"0 0 30 0 0 ML\" #Left-justify so montage titles are uniform.\n");
      fprintf(new_fp,"  title_x=$(bc <<< \"scale=5; $blurb_x-0.1\")\n");
      fprintf(new_fp,"else\n");
      fprintf(new_fp,"  title_x=$(bc <<< \"scale=5; $map_x+$map_width/2\")\n");
      fprintf(new_fp,"fi\n");
      fprintf(new_fp,"title_pos=\" -Xa${title_x}c -Ya${title_y}c \"\n");
      fprintf(new_fp,"echo $title_format $title | $gmt_prefix pstext -N $title_pos $misc_range $middle >> $plot_base.ps\n");
      
      //Draw color scale with units printed above.
      write_gmt_colorscale(new_fp,results,0);//0 = next to map for GMT PDF output.

      //Print blurb about data range, or masked amplitudes for phase plots.
      //Don't print blurb if disabled or if this is an RGB map.
      if(results.rgb.size()==3 and results.latlon.outputs.size()==3) plot_options.blurb_disabled = 1;
      if(plot_options.blurb_disabled) fprintf(new_fp,"blurb_contents=\"\"\n");
      fprintf(new_fp,"echo $blurb_format $blurb_contents | $gmt_prefix pstext -N $blurb_pos $misc_range $middle >> $plot_base.ps\n");

      //Only print error bars if they're the right length, and this one is > 0.0.
      //Need to have separate copies for unstructured and latlon maps because
      //in each case the "right length" is defined by a different multivector.
      int blurb2_written=0;
      if(results.latlon.outputs.size() == results.error_bars.size() and !results.latlon.outputs.empty()){
        if(results.error_bars.at(i) > 0.0){
          blurb2_written=1;
          fprintf(new_fp,"blurb2_contents=\"Error bar: %.1f $scale_units\"\n",results.error_bars.at(i));
        }
      }
      if(!blurb2_written) fprintf(new_fp,"blurb2_contents=\"\" #Error bar: N/A $scale_units\"\n");
      fprintf(new_fp,"echo $blurb_format $blurb2_contents | $gmt_prefix pstext -N $blurb2_pos $misc_range $end >> $plot_base.ps\n");
    }
    else cout<<"!!!!WARNING!!!!!! results.options.output_choice "<<results.options.output_choice<<" isn't recognized."<<endl;
  }
  else{
    if(results.options.output_choice == 5){
      //Plot ONLY the data in cylindrical equidistant projection.
      //Skip title, but draw coastlines b/c Google's coastline is a very thin, hard-to-see yellow line... AND it's inaccurate in Antarctica.
      if(results.rgb.size()==3 and results.latlon.outputs.size()==3) fprintf(new_fp,"$gmt_prefix grdimage ../red.nc ../green.nc ../blue.nc $resolution $range $projection $map_pos $start > $plot_base.ps\n");
      else fprintf(new_fp,"$gmt_prefix grdimage $data_name $resolution $range $projection $map_pos -Cmap.cpt $start > $plot_base.ps\n");
      write_gmt_coastlines(new_fp);
      //If marker_lats multivectors have the right size, plot markers in the color "sandy brown".
      //Plot markers before mascons because mascons are smaller than markers.
      if(results.marker_lats.size() == results.latlon.outputs.size() and results.marker_lons.size() == results.latlon.outputs.size() and !results.latlon.outputs.empty()){
        if(!results.marker_lats.at(i).empty() and !results.marker_lons.at(i).empty()){
          fprintf(new_fp,"$gmt_prefix psxy -N $data_name -bcmarker_lons/marker_lats -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n");
        }
      }
      //Plot mascon centers in the color "saddle brown".
      //DISABLED: mascon centers are distorted near the poles in Google Earth.
      //fprintf(new_fp,"$gmt_prefix psxy $data_name -bcmascon_lons/mascon_lats -Sc0.01c -G139/69/19 $range $projection $map_pos $middle >> $plot_base.ps\n");
      fprintf(new_fp,"#So that previous statements can all use $middle rather than requiring\n");
      fprintf(new_fp,"#if-then statements to separate cases requiring $end.\n");
      fprintf(new_fp,"echo $blurb_format | $gmt_prefix pstext -N $blurb_pos $misc_range $end >> $plot_base.ps\n");
      //Print blurb about data range, or masked amplitudes for phase plots.
      //Because this is a KMZ plot, need to write it by itself.
      if(plot_options.blurb_disabled) fprintf(new_fp,"#");
      fprintf(new_fp,"echo $blurb_format $blurb_contents | $gmt_prefix pstext -N $blurb_pos $misc_range > blurb_$plot_base.ps\n");
      //Print title.
      //Because this is a KMZ plot, need to write it by itself.
      fprintf(new_fp,"title=\"%s\"\n",results.titles.at(i).c_str());
      fprintf(new_fp,"echo $blurb_format $title | $gmt_prefix pstext -N $blurb_pos $misc_range > title_$plot_base.ps\n");
    }
    else cout<<"!!!!WARNING!!!!!! results.options.output_choice "<<results.options.output_choice<<" isn't recognized."<<endl;
  }
}

void write_gmt_plot_data(results_s &results,
                         FILE *new_fp,
                         long long i){
  /**********************************************************************
  Purpose: This function writes the GMT data plotting command for 2D plots.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          new_fp - file pointer to current gmt script file
          i - (output) index of parameter being plotted.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long j;

  //Search for min/max, store in results.xy.x/ymin,x/ymax.
  minmax_output(results,i,0);//Map index, quiet.

  //GMT fails if max==min.
  double margin = 0.1;
  if(results.xy.xmin == results.xy.xmax){
    cout<<"NOTE: Both xmin and xmax = "<<results.xy.xmin<<" so they're each being shifted by "<<margin*100<<"%."<<endl;
    if(fabs(results.xy.xmin) > tolerance){
      results.xy.xmin -= margin*fabs(results.xy.xmin);
      results.xy.xmax += margin*fabs(results.xy.xmax);
    }
    else{
      results.xy.xmin -= margin;
      results.xy.xmax += margin;
    }
  }
  if(results.xy.ymin == results.xy.ymax){
    cout<<"NOTE: Both ymin and ymax = "<<results.xy.ymin<<" so they're each being shifted by "<<margin*100<<"%."<<endl;
    if(fabs(results.xy.ymin) > tolerance){
      results.xy.ymin -= margin*fabs(results.xy.ymin);
      results.xy.ymax += margin*fabs(results.xy.ymax);
    }
    else{
      results.xy.ymin -= margin;
      results.xy.ymax += margin;
    }
  }
  
  //Specify how many tick marks are on the x/y axes.
  long long num_ticks = 10;
  double x_tick_spacing = (results.xy.xmax-results.xy.xmin)/(double)num_ticks;
  double y_tick_spacing = (results.xy.ymax-results.xy.ymin)/(double)num_ticks;

  fprintf(new_fp,"line_thk=\"5\"\n");
  fprintf(new_fp,"sym_size=\"0.3c\"\n");
  fprintf(new_fp,"range=\" -R%.10f/%.10f/%.10f/%.10f \"\n",results.xy.xmin,results.xy.xmax,results.xy.ymin,results.xy.ymax);
  fprintf(new_fp,"plot_pos=\"-Xa4.5c -Ya3c\"\n");
  fprintf(new_fp,"legend_pos=\"-Dx5c/8.5c/13c/9.0c/BL\"\n");
  string x_format("24c");
  string y_format("15c");
  if(results.xy.log_x_axis.size() == results.xy.x_values.size()){
    if(results.xy.log_x_axis.at(i) == 1) x_format.append("l");
  }
  if(results.xy.log_y_axis.size() == results.xy.y_values.size()){
    if(results.xy.log_y_axis.at(i) == 1) y_format.append("l");
  }
  fprintf(new_fp,"projection=\" -JX%s/%s \"\n",x_format.c_str(),y_format.c_str());
  fprintf(new_fp,"main_title=\"%s\"\n",results.xy.titles.at(i).c_str());
  //A more automated way to set sigfigs would be to take logarithm of spacing,
  //then use setprecision() and c++ output routines to write digits using
  //the logarithm.
  fprintf(new_fp,"x_axis=\"a%.8f\"\n",x_tick_spacing);
  fprintf(new_fp,"y_axis=\"a%.8f\"\n",y_tick_spacing);
  fprintf(new_fp,"x_units=\"%s\"\n",results.xy.x_units.at(i).c_str());
  fprintf(new_fp,"y_units=\"%s\"\n",results.xy.y_units.at(i).c_str());
  fprintf(new_fp,"title_pos=\" -Xa0c -Ya19c \"\n");
  fprintf(new_fp,"title_format=\"0 0 30 0 0 ML\"\n");
  fprintf(new_fp,"montagetitle=\"\"\n");
  fprintf(new_fp,"if [ $montage != 0 ]\n");
  fprintf(new_fp,"then\n");
  fprintf(new_fp,"  montagetitle=${prefixes[$index]}\" \"$main_title\n");
  fprintf(new_fp,"  main_title=\"\"\n");
  fprintf(new_fp,"fi\n");
  fprintf(new_fp,"#This puts 4 digits on y-axis, 0 digits on x-axis.\n");
  fprintf(new_fp,"#$gmt_prefix gmtset D_FORMAT=%%.4f\n");
  fprintf(new_fp,"#$gmt_prefix psbasemap $plot_pos $range $projection -B/$y_axis:\"$y_units\"::.\"$main_title\":Wesn $start > $plot_base.ps\n");
  fprintf(new_fp,"#$gmt_prefix gmtset D_FORMAT=%%.0f\n");
  fprintf(new_fp,"#$gmt_prefix psbasemap $plot_pos $range $projection -B$x_axis:\"$x_units\":weSn $middle >> $plot_base.ps\n");
  fprintf(new_fp,"$gmt_prefix psbasemap $plot_pos $range $projection -B$x_axis:\"$x_units\":/$y_axis:\"$y_units\"::.\"$main_title\":WeSn $start > $plot_base.ps\n");
  fprintf(new_fp,"echo $title_format $montagetitle | $gmt_prefix pstext -N $title_pos -R0/21.59/0/27.94 -Jx1c $middle >> $plot_base.ps\n");
  //Define variable with all NetCDF variable names for different lines in this plot.
  fprintf(new_fp,"all_lines=\"");
  for(j=0;j<(long long)results.xy.x_values.at(i).size();j++){
    fprintf(new_fp,"x_values_%04lld/",j+1);
    fprintf(new_fp,"y_values_%04lld ",j+1);
  }
  fprintf(new_fp,"\"\n");//Finish list of NetCDF variables and endline.
  //If legends are present for all lines, define them.
  int legends = 0;
  fprintf(new_fp,"all_legends=(");
  if((long long)results.xy.legends.size() > i){
    if(results.xy.legends.at(i).size() == results.xy.x_values.at(i).size()){
      legends = 1;
      for(j=0;j<(long long)results.xy.x_values.at(i).size();j++){
        fprintf(new_fp,"\"%s\" ",results.xy.legends.at(i).at(j).c_str());
      }
    }
  }
  if(legends==0) for(j=0;j<(long long)results.xy.x_values.at(i).size();j++) fprintf(new_fp,"\"\" ");
  fprintf(new_fp,")\n");//Finish array of legends and endline.

  fprintf(new_fp,"all_colors=(\"black\" \"red\" \"green\" \"blue\")\n");
  fprintf(new_fp,"all_symbols=(\"c\" \"d\" \"s\" \"t\" \"+\" \"x\")\n");
  fprintf(new_fp,"num_colors=${#all_colors[*]}\n");
  fprintf(new_fp,"num_symbols=${#all_symbols[*]}\n");
  fprintf(new_fp,"i=0\n");
  fprintf(new_fp,"echo \"#This file contains information used by pslegend.\" > legend.txt\n");
  fprintf(new_fp,"for current_line in $all_lines\n");
  fprintf(new_fp,"do\n");
  fprintf(new_fp,"  let \"col=$i %% $num_colors\"\n");
  fprintf(new_fp,"  let \"sym=$i %% $num_symbols\"\n");
  fprintf(new_fp,"  #Plot lines.\n");
  fprintf(new_fp,"  $gmt_prefix psxy -N $data_name -bc$current_line -W$line_thk/${all_colors[$col]} $range $projection $plot_pos $middle >> $plot_base.ps\n");
  fprintf(new_fp,"  #Plot symbols.\n");
  fprintf(new_fp,"  $gmt_prefix psxy -N $data_name -bc$current_line -S${all_symbols[$sym]}$sym_size -G${all_colors[$col]} $range $projection $plot_pos $middle >> $plot_base.ps\n");
  fprintf(new_fp,"  #Record legend information.\n");
  fprintf(new_fp,"  echo \"S 0.1c ${all_symbols[$sym]} $sym_size ${all_colors[$col]} 0.25p 0.8c ${all_legends[$i]}\" >> legend.txt\n");
  fprintf(new_fp,"  let \"i++\"\n");
  fprintf(new_fp,"done\n");
  if(legends==0) fprintf(new_fp,"#");
  fprintf(new_fp,"$gmt_prefix pslegend legend.txt $range $projection $legend_pos $middle >> $plot_base.ps\n");
  fprintf(new_fp,"#Uncomment to put a red arrow w black border pointing to \"x y angle size\".\n");
  fprintf(new_fp,"#echo 18.0 505.329954 45 1 | $gmt_prefix psxy -N -Svh -Wblack -Gred $range $projection $plot_pos $middle >> $plot_base.ps\n");
  //Loop through series of X,Y labels, but only if they're recorded for this plot.
  if((long long)results.xy.x_label_values.size() > i){//Because i might be too large!
    if(!results.xy.x_label_values.at(i).empty()){//Only bother if there ARE labels.
      fprintf(new_fp,"offset=\"-D0.3cv\"\n");
      fprintf(new_fp,"label_format=\"15 0 1 2\"\n");
      for(j=0;j<(long long)results.xy.x_label_values.at(i).size();j++){
        fprintf(new_fp,"echo %.10f %.10f $label_format %s | $gmt_prefix pstext -N $plot_pos $offset $range $projection $middle >> $plot_base.ps\n",results.xy.x_label_values.at(i).at(j),results.xy.y_label_values.at(i).at(j),results.xy.labels.at(i).at(j).c_str());
      }
    }
  }
  if(results.rgb.size() == 3){
    fprintf(new_fp,"#Plot square with rgb color.\n");
    fprintf(new_fp,"echo 21 12 | $gmt_prefix psxy -N -Ss4.0c -G%d/%d/%d -W0.1c,black -R0/21.59/0/27.94 -Jx1c $plot_pos $middle >> $plot_base.ps\n",(int)(results.rgb.at(0)*255.0),(int)(results.rgb.at(1)*255.0),(int)(results.rgb.at(2)*255.0));
  }
  fprintf(new_fp,"#So that previous statements can all use $middle rather than requiring\n");
  fprintf(new_fp,"#if-then statements to separate cases requiring $end.\n");
  fprintf(new_fp,"echo 0 90 15 0 1 ML | $gmt_prefix pstext -N -Xa-10.5c -Ya-7c $range $projection $end >> $plot_base.ps\n");
}

void write_gmt_coastlines(FILE *new_fp){
  /**********************************************************************
  Purpose: This function writes the GMT coastlines command(s).
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          new_fp - file pointer to current gmt script file
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  fprintf(new_fp,"$gmt_prefix pscoast -W$pscoast_thk/$coast_color $pscoast_res $range $projection $map_pos $middle >> $plot_base.ps\n");
  fprintf(new_fp,"if [ $coastlines == 2 ]\n");
  fprintf(new_fp,"then\n");
  fprintf(new_fp,"  $gmt_prefix psxy -N $coast_file -: -Sc$coast_thk -W$coast_thk/$coast_color $range $projection $map_pos $middle >> $plot_base.ps\n");
  fprintf(new_fp,"fi\n");
}

void write_gmt_colorscale(FILE *new_fp,
                          results_s &results,
                          int kml_output){
  /**********************************************************************
  Purpose: This function writes the GMT color scale command.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          new_fp - file pointer to current gmt script file
          kml_output - 0 for GMT scale next to plot, 1 for KMZ by itself.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  results_s rgb;//Just to init choice and maxes.

  if(results.rgb.size() != 3) fprintf(new_fp,"cpt_name=\"-Cmap.cpt \"\n");
  else fprintf(new_fp,"cpt_name=\"-C../../rgb00001.cpt \"\n");
  if(results.rgb.size() != 3 or rgb.rgb_choice < 2){
    //If this is intended for KMZ output, don't overlay the scale.
    if(kml_output){
      //KMZ version writes to a different file.
      fprintf(new_fp,"$gmt_prefix psscale $cpt_name -L $scale_format $overflow $scale_pos -A $start > scale_$plot_base.ps\n");
      fprintf(new_fp,"#Print units manually, otherwise they're too close to numbers on scale.\n");
      fprintf(new_fp,"echo $units_format $scale_units | $gmt_prefix pstext -N $units_pos $misc_range $end >> scale_$plot_base.ps\n");
    }
    else{
      fprintf(new_fp,"$gmt_prefix psscale $cpt_name -L $scale_format $overflow $scale_pos -A $middle >> $plot_base.ps\n");
      fprintf(new_fp,"#Print units manually, otherwise they're too close to numbers on scale.\n");
      fprintf(new_fp,"echo $units_format $scale_units | $gmt_prefix pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n");
    }
  }
  else{
    fprintf(new_fp,"numwidths=%lld\n",rgb.max_widths);
    fprintf(new_fp,"$gmt_prefix gmtset TICK_LENGTH 0.3c\n");
    fprintf(new_fp,"scale_width=$(bc <<< \"scale=5; $scale_width / $numwidths\")\n");
    fprintf(new_fp,"for (( j=1; j <= $numwidths; j++ ))\n");
    fprintf(new_fp,"do\n");
    if(kml_output){
      //KMZ version writes to a different file.
      fprintf(new_fp,"  $gmt_prefix psscale $cpt_name -L $scale_format $overflow $scale_pos -A $start > scale_$plot_base.ps\n");
      fprintf(new_fp,"  #Print units manually, otherwise they're too close to numbers on scale.\n");
      fprintf(new_fp,"  echo $units_format $scale_units | $gmt_prefix pstext -N $units_pos $misc_range $end >> scale_$plot_base.ps\n");
    }
    else{
      fprintf(new_fp,"  j_string=$(printf '%%05d' $j)\n");
      fprintf(new_fp,"  cpt_name=\"-C../rgb$j_string.cpt\"\n");
      fprintf(new_fp,"  $gmt_prefix psscale $cpt_name -L $scale_format $overflow $scale_pos -S -A $middle >> $plot_base.ps\n");
      fprintf(new_fp,"  #Print units manually, otherwise they're too close to numbers on scale.\n");
      fprintf(new_fp,"  echo $units_format $scale_units | $gmt_prefix pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n");
      fprintf(new_fp,"  #Move next scale to the right and get rid of the tick marks.\n");
      fprintf(new_fp,"  scale_x=$(bc <<< \"scale=5; $scale_x+$scale_width\")\n");
      fprintf(new_fp,"  scale_pos=\" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c \"\n");
      fprintf(new_fp,"  $gmt_prefix gmtset TICK_LENGTH 0.0\n");
    }
    fprintf(new_fp,"done\n");
  }
}

void write_rgb_colorscale(results_s &results, results_s &cie, plot_options_s &plot_options, int verbose){
  /**********************************************************************
  Purpose: This function writes the RGB colorscale to outputfolder.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long i,j,l;
  FILE *new_fp;
  char s[max_length];
  string new_file;
  results_s old,copy,all;

  //rgb_choice and maxes initialized in definitions.hpp instead of here.

  switch(results.rgb_choice){
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    case 0:{
      //Create CPT file in outputfolder, which is NOT netcdf_output bc that doesn't exist yet.
      new_file = plot_options.outputfolder;
      new_file.append("rgb00001.cpt");
      new_fp = fopen(new_file.c_str(),"w");
      //Warn if file doesn't open correctly.
      if(!new_fp) cout<<"The rgb00001.cpt file couldn't be created."<<endl;
      
      fprintf(new_fp,"# COLOR_MODEL = RGB\n");
      cout<<"Writing RGB colorscale to disk."<<endl;
      //Loop through periods for colorscale divisions.
      for(i=0;i<(long long)results.xy.x_values.at(0).at(0).size();i++){
        //Calculate color for each period.
        copy = results;
        create_synthetic_plot(copy,20,0,0,0,results.xy.x_values.at(0).at(0).at(i),0.2);//choice,#pts,x0,delta,param1,2.
        match_x_axes(copy,cie);
        interpolate(copy,cie);
        spectrum2xyz(copy,cie,verbose);
        xyz2rgb(copy);
        rescale_single_rgb(copy,verbose);
        gamma_correct_single_rgb(copy,verbose);
        for(l=0;l<3;l++) copy.rgb.at(l) *= 255.0;
        if(i>0) fprintf(new_fp,"%12.6e %3d %3d %3d %12.6e %3d %3d %3d\n",results.xy.x_values.at(0).at(0).at(i-1),(int)old.rgb.at(0),(int)old.rgb.at(1),(int)old.rgb.at(2),results.xy.x_values.at(0).at(0).at(i),(int)copy.rgb.at(0),(int)copy.rgb.at(1),(int)copy.rgb.at(2));
        old = copy;
      }
      fclose(new_fp);
      break;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    case 1:{
      //Create CPT file in outputfolder, which is NOT netcdf_output bc that doesn't exist yet.
      new_file = plot_options.outputfolder;
      new_file.append("rgb00001.cpt");
      new_fp = fopen(new_file.c_str(),"w");
      //Warn if file doesn't open correctly.
      if(!new_fp) cout<<"The rgb00001.cpt file couldn't be created."<<endl;
      
      fprintf(new_fp,"# COLOR_MODEL = RGB\n");
      cout<<"Writing RGB colorscale to disk."<<endl;
      //Shortest (first) periods are spaced closest together. Smaller half-widths aren't needed.
      double hw = results.xy.x_values.at(0).at(0).at(1) - results.xy.x_values.at(0).at(0).at(0);
      hw*=2;/*Otherwise 2-5 year period colorscales fail like this:
      (X,Y,Z) = (5.3566964e+09, -1.7371269e+07, 2.6757576e+10)
      Initial (r,g,b)= (-3.4607248e+07, -1.6634301e-05, -1.3869090e+08)
      Rescaling by -1.6634300790e-05
      Linear (r,g,b) = (2.0804751e+12, 1.0000000e+00, 8.3376450e+12)
      Nonlinear (r,g,b) = (3.8391046e+05, 1.0000000e+00, 7.1701195e+05)
      (X,Y,Z) = (1.1417463e+10, 2.9056321e+09, 5.4360743e+10)
      Initial (r,g,b)= (4.0777485e+09, 1.8195913e-03, 2.8237678e+10)
      Rescaling by 2.8237678453e+10
      Linear (r,g,b) = (1.4440807e-01, 6.4438418e-14, 1.0000000e+00)
      Nonlinear (r,g,b) = (3.6097219e-01, 2.9041191e-13, 1.0000000e+00)
      2.000000e+00 97897166 255 182838046 2.250137e+00  92   0 255 // */
      all.options.output_choice = 5;
      all.latlon.lat.resize(results.xy.x_values.at(0).at(0).size());
      all.latlon.lon.resize(1);
      j=0;//Because there's only 1 lon.
      init_latlon(all,3);
      //Loop through periods for colorscale divisions.
      for(i=0;i<(long long)results.xy.x_values.at(0).at(0).size();i++){
        //Calculate color for each period.
        copy = results;
        create_synthetic_plot(copy,20,0,0,0,results.xy.x_values.at(0).at(0).at(i),hw);//choice,#pts,x0,delta,param1,2.
        match_x_axes(copy,cie);
        interpolate(copy,cie);
        spectrum2xyz(copy,cie,verbose);
        xyz2rgb(copy);
        for(l=0;l<3;l++) all.latlon.outputs.at(l).at(i).at(j) = copy.rgb.at(l);
      }
      spectra_end(all,2000,1);//mod_choice,norm_rgb
      for(i=1;i<(long long)results.xy.x_values.at(0).at(0).size();i++){
        fprintf(new_fp,"%12.6e %3d %3d %3d %12.6e %3d %3d %3d\n",results.xy.x_values.at(0).at(0).at(i-1),(int)all.latlon.outputs.at(0).at(i-1).at(j),(int)all.latlon.outputs.at(1).at(i-1).at(j),(int)all.latlon.outputs.at(2).at(i-1).at(j),results.xy.x_values.at(0).at(0).at(i),(int)all.latlon.outputs.at(0).at(i).at(j),(int)all.latlon.outputs.at(1).at(i).at(j),(int)all.latlon.outputs.at(2).at(i).at(j));
      }
      fclose(new_fp);
      break;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    case 2:{
      cout<<"Writing RGB colorscale to disk."<<endl;
      long long label_interval;
      //Shortest (first) periods are spaced closest together. Smaller half-widths aren't needed.
      double hw = results.xy.x_values.at(0).at(0).at(1) - results.xy.x_values.at(0).at(0).at(0);
      //hw*=2;//Otherwise 2-5 year period colorscales fail as documented in case 1.
      //hw*=4;//Otherwise first colorbar is very dim and the blue side has a lot of black.
      double hw_end = 0.35*(results.xy.x_values.at(0).at(0).back() - results.xy.x_values.at(0).at(0).at(0));
      if(hw_end <= hw) hw_end = (results.xy.x_values.at(0).at(0).back() - results.xy.x_values.at(0).at(0).at(0));
      double hw_step = (hw_end-hw)/(results.max_widths-1);
      if((long long)results.xy.x_values.at(0).at(0).size() > results.max_labels) label_interval = (long long)((double)results.xy.x_values.at(0).at(0).size() / results.max_labels);
      else label_interval = 1;
      all.options.output_choice = 5;
      all.latlon.lat.resize(results.xy.x_values.at(0).at(0).size());
      all.latlon.lon.resize(results.max_widths);
      init_latlon(all,3);
      
      for(j=0;j<results.max_widths;j++){
        //Loop through periods for colorscale divisions.
        for(i=0;i<(long long)results.xy.x_values.at(0).at(0).size();i++){
          //Calculate color for each period.
          copy = results;
          create_synthetic_plot(copy,20,0,0,0,results.xy.x_values.at(0).at(0).at(i),hw);//choice,#pts,x0,delta,param1,2.
          match_x_axes(copy,cie);
          interpolate(copy,cie);
          spectrum2xyz(copy,cie,verbose);
          xyz2rgb(copy);
          for(l=0;l<3;l++) all.latlon.outputs.at(l).at(i).at(j) = copy.rgb.at(l);
        }
        hw+=hw_step;
      }
      spectra_end(all,2000,1);//mod_choice,norm_rgb
      for(j=0;j<results.max_widths;j++){
        //Create CPT file in outputfolder, which is NOT netcdf_output bc that doesn't exist yet.
        new_file = plot_options.outputfolder;
        sprintf(s,"rgb%05lld.cpt",j+1);
        new_file.append(s);
        new_fp = fopen(new_file.c_str(),"w");
        //Warn if file doesn't open correctly.
        if(!new_fp) cout<<"The rgb.cpt file couldn't be created."<<endl;
        fprintf(new_fp,"# COLOR_MODEL = RGB\n");
        for(i=1;i<(long long)results.xy.x_values.at(0).at(0).size();i++){
          fprintf(new_fp,"%12.6e %3d %3d %3d %12.6e %3d %3d %3d",results.xy.x_values.at(0).at(0).at(i-1),(int)all.latlon.outputs.at(0).at(i-1).at(j),(int)all.latlon.outputs.at(1).at(i-1).at(j),(int)all.latlon.outputs.at(2).at(i-1).at(j),results.xy.x_values.at(0).at(0).at(i),(int)all.latlon.outputs.at(0).at(i).at(j),(int)all.latlon.outputs.at(1).at(i).at(j),(int)all.latlon.outputs.at(2).at(i).at(j));
          if(j>0) fprintf(new_fp," ;\n");
          else if((i-1) % label_interval==0) fprintf(new_fp," ; %.1f\n",results.xy.x_values.at(0).at(0).at(i-1));
          else fprintf(new_fp," ;\n");
        }
        fclose(new_fp);
      }
      break;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    default:cout<<"!!!WARNING!!! Didn't recognize results.rgb_choice "<<results.rgb_choice<<endl;
  }//End of switch-case-break.
}


void write_gmt_scripts(plot_options_s &plot_options,
                       grid_s &grid,
                       results_s &results){
  /**********************************************************************
  Purpose: This function writes the GMT scripts.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long i;
  char s[max_length];
  FILE *new_fp, *flip_fp, *trim_fp, *extra_fp;
  string new_file, flip_file, trim_file, extra_file, temp_string;

  string output_base;

  //Create GMT script in outputfolder, which should be netcdf_output.
  new_file = plot_options.outputfolder;
  new_file.append("create_plots.sh");
  new_fp = fopen(new_file.c_str(),"w");
  //Warn if file doesn't open correctly.
  if(!new_fp) cout<<"The create_plots.sh GMT script couldn't be created."<<endl;
  
  //Every bash script needs this to be the first line.
  fprintf(new_fp,"#!/bin/bash\n");
  fprintf(new_fp,"#set -x #Uncomment to echo these commands.\n");

  //Create flip_backgrounds script in outputfolder, which should be netcdf_output.
  flip_file = plot_options.outputfolder;
  flip_file.append("flip_backgrounds.sh");
  flip_fp = fopen(flip_file.c_str(),"w");
  //Warn if file doesn't open correctly.
  if(!flip_fp) cout<<"The flip_backgrounds.sh script couldn't be created."<<endl;
  //Every bash script needs this to be the first line.
  fprintf(flip_fp,"#!/bin/bash\n");
  fprintf(flip_fp,"set -x\n");

  //Create trim script in outputfolder, which should be netcdf_output.
  trim_file = plot_options.outputfolder;
  trim_file.append("trim.sh");
  trim_fp = fopen(trim_file.c_str(),"w");
  //Warn if file doesn't open correctly.
  if(!trim_fp) cout<<"The trim.sh script couldn't be created."<<endl;
  //Every bash script needs this to be the first line.
  fprintf(trim_fp,"#!/bin/bash\n");
  fprintf(trim_fp,"set -x\n");

  switch(results.options.output_choice){
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    case 1:{//Unstructured map output format.
      output_base = "map_parameter";
      for(i=0;i<(long long)plot_options.just_the_filenames.size();i++){
        //If this is the first file...
        if(i==0){
          write_gmt_defs(new_fp);
          fprintf(new_fp,"color_scheme=%d #1/2=white/black background.\n",plot_options.color_scheme);
          fprintf(new_fp,"montage=%d #1=left-justify titles, add (a),(b), run montage.sh.\n",plot_options.montage);
          fprintf(new_fp,"prefixes=('(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' '(k)' '(l)' '(m)' '(n)' '(o)' '(p)' '(q)' '(r)' '(s)' '(t)' '(u)' '(v)' '(w)' '(x)' '(y)' '(z)')\n");
          fprintf(new_fp,"index=-1 #Increments on each map, accesses prefixes above for montage.\n");
          fprintf(new_fp,"png_options=\" -P -Tg \" #PDF default: -E720, else 300 dpi.\n");
          fprintf(new_fp,"if [ $montage != 0 ]\n");
          fprintf(new_fp,"then\n");
          fprintf(new_fp,"  png_options=\" -A\"$png_options\n");
          fprintf(new_fp,"fi\n");
          fprintf(new_fp,"#Force red foreground color because maps with points are white.\n");
          fprintf(new_fp,"$gmt_prefix gmtset COLOR_BACKGROUND=2/2/2 COLOR_FOREGROUND=red\n");
          fprintf(new_fp,"digits=%d\n",plot_options.scale_digits);
          fprintf(new_fp,"$gmt_prefix gmtset D_FORMAT=%%.${digits}f\n");
        }
        //Define the filenames for PostScript output.
        sprintf(s,"%s_%04lld",output_base.c_str(),i+1);
        fprintf(new_fp,"#######################################################\n");
        fprintf(new_fp,"data_name=\"%s\"\n",plot_options.just_the_filenames.at(i).c_str());
        fprintf(new_fp,"plot_base=\"%s\"\n",s);
        fprintf(new_fp,"let index=$index+1\n");
        fprintf(new_fp,"#######################################################\n");
        //Record the filename bases in the flip_backgrounds.sh and trim.sh scripts.
        if(i==0){
          if(plot_options.just_the_filenames.size()==1){
            fprintf(flip_fp,"all_bases=\"%s\"\n",s);//First and last base starts and ends the string list.
            fprintf(trim_fp,"all_bases=\"%s\"\n",s);//First and last base starts and ends the string list.
          }
          else{
            fprintf(flip_fp,"all_bases=\"%s\n",s);//First base starts the string list.
            fprintf(trim_fp,"all_bases=\"%s\n",s);//First base starts the string list.
          }
        }
        else if(i<(long long)plot_options.just_the_filenames.size()-1){
          fprintf(flip_fp,"%s\n",s);//Bases that aren't first or last are unadorned.
          fprintf(trim_fp,"%s\n",s);//Bases that aren't first or last are unadorned.
        }
        else{
          fprintf(flip_fp,"%s\"\n",s);//Last base ends with a " to terminate string list.
          fprintf(trim_fp,"%s\"\n",s);//Last base ends with a " to terminate string list.
        }

        //Plot data, with title on top and coastlines.
        write_gmt_map_data(results,grid,plot_options,new_fp,results.titles.at(i),i);

        fprintf(new_fp,"ps2raster $png_options $plot_base.ps #Convert PS to PNG format.\n");
        fprintf(new_fp,"#ps2raster -P -Tf $plot_base.ps #Convert PS to PDF, if uncommented.\n");

        //Delete PS file because it's twice as large (more at -E2000 and 0.25x0.25 global- 2.7MB PDF, 165MB PS!) as the PDF or PNG.
        fprintf(new_fp,"rm -f $plot_base.ps\n");
        
        //Move cpt file to backup version so it won't be automatically used
        //if this script is executed again... just in case write_gmt_cpt fails:
        //it's good for that failure to be obvious.
        if(results.rgb.size() != 3) fprintf(new_fp,"mv map.cpt Zbackup_cpt_$plot_base.cpt\n");
        
        //When this plot is finished, store its data_name in previous_data_name
        //so that the next plot (if it's a phase plot with amplitude masking) can
        //access that data.
        fprintf(new_fp,"previous_data_name=$data_name\n");
      }
      break;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    case 2:{//2D plot output format.
      
      output_base = "plot_parameter";
      for(i=0;i<(long long)plot_options.just_the_filenames.size();i++){
        //If this is the first file...
        if(i==0){
          write_gmt_defs(new_fp);
          fprintf(new_fp,"color_scheme=%d\n",plot_options.color_scheme);
          fprintf(new_fp,"montage=%d #1=left-justify titles, add (a),(b), run montage.sh.\n",plot_options.montage);
          fprintf(new_fp,"prefixes=('(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' '(k)' '(l)' '(m)' '(n)' '(o)' '(p)' '(q)' '(r)' '(s)' '(t)' '(u)' '(v)' '(w)' '(x)' '(y)' '(z)')\n");
          fprintf(new_fp,"index=-1 #Increments on each map, accesses prefixes above for montage.\n");
          fprintf(new_fp,"png_options=\" -P -Tg \" #PDF default: -E720, else 300 dpi.\n");
          fprintf(new_fp,"if [ $montage != 0 ]\n");
          fprintf(new_fp,"then\n");
          fprintf(new_fp,"  png_options=\" -A\"$png_options\n");
          fprintf(new_fp,"fi\n");
          fprintf(new_fp,"#Force off-white(dark gray) fore(back)ground color because\n");
          fprintf(new_fp,"#flip_backgrounds.sh can change the plots' text from\n");
          fprintf(new_fp,"#black to white, and their backgrounds from white to black.\n");
          fprintf(new_fp,"$gmt_prefix gmtset COLOR_BACKGROUND=2/2/2 COLOR_FOREGROUND=253/253/253\n");
          fprintf(new_fp,"digits=%d\n",plot_options.scale_digits);
          fprintf(new_fp,"$gmt_prefix gmtset D_FORMAT=%%.${digits}f\n");
        }
        //Define the filenames for PostScript output.
        sprintf(s,"%s_%04lld",output_base.c_str(),i+1);
        fprintf(new_fp,"#######################################################\n");
        fprintf(new_fp,"data_name=\"%s\"\n",plot_options.just_the_filenames.at(i).c_str());
        fprintf(new_fp,"plot_base=\"%s\"\n",s);
        fprintf(new_fp,"let index=$index+1\n");
        fprintf(new_fp,"#######################################################\n");
        //Record the filename bases in the flip_backgrounds.sh and trim.sh scripts.
        if(i==0){
          if(plot_options.just_the_filenames.size()==1){
            fprintf(flip_fp,"all_bases=\"%s\"\n",s);//First and last base starts and ends the string list.
            fprintf(trim_fp,"all_bases=\"%s\"\n",s);//First and last base starts and ends the string list.
          }
          else{
            fprintf(flip_fp,"all_bases=\"%s\n",s);//First base starts the string list.
            fprintf(trim_fp,"all_bases=\"%s\n",s);//First base starts the string list.
          }
        }
        else if(i<(long long)plot_options.just_the_filenames.size()-1){
          fprintf(flip_fp,"%s\n",s);//Bases that aren't first or last are unadorned.
          fprintf(trim_fp,"%s\n",s);//Bases that aren't first or last are unadorned.
        }
        else{
          fprintf(flip_fp,"%s\"\n",s);//Last base ends with a " to terminate string list.
          fprintf(trim_fp,"%s\"\n",s);//Last base ends with a " to terminate string list.
        }

        //Plot data.
        write_gmt_plot_data(results,new_fp,i);

        fprintf(new_fp,"ps2raster $png_options $plot_base.ps #Convert PS to PNG format.\n");
        fprintf(new_fp,"#ps2raster -P -Tf $plot_base.ps #Convert PS to PDF, if uncommented.\n");

        //Delete PS file because it's twice as large (more at -E2000 and 0.25x0.25 global- 2.7MB PDF, 165MB PS!) as the PDF or PNG.
        fprintf(new_fp,"rm -f $plot_base.ps\n");
      }
      break;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    case 3:{//Maps of separated points- no contours as in types 0,1.
      cout<<"!!!WARNING!!! Type "<<results.options.output_choice<<" doesn't have GMT script routines written yet!"<<endl;
      break;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    case 4:{//Maps of separated points with labels.

      output_base = "map_number";
      for(i=0;i<(long long)plot_options.just_the_filenames.size();i++){
        //If this is the first file...
        if(i==0){
          write_gmt_defs(new_fp);
          fprintf(new_fp,"color_scheme=%d #1/2=white/black background\n",plot_options.color_scheme);
          fprintf(new_fp,"montage=%d #1=left-justify titles, add (a),(b), run montage.sh.\n",plot_options.montage);
          fprintf(new_fp,"prefixes=('(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' '(k)' '(l)' '(m)' '(n)' '(o)' '(p)' '(q)' '(r)' '(s)' '(t)' '(u)' '(v)' '(w)' '(x)' '(y)' '(z)')\n");
          fprintf(new_fp,"index=-1 #Increments on each map, accesses prefixes above for montage.\n");
          fprintf(new_fp,"png_options=\" -P -Tg \" #PDF default: -E720, else 300 dpi.\n");
          fprintf(new_fp,"if [ $montage != 0 ]\n");
          fprintf(new_fp,"then\n");
          fprintf(new_fp,"  png_options=\" -A\"$png_options\n");
          fprintf(new_fp,"fi\n");
          fprintf(new_fp,"#Force off-white(dark gray) fore(back)ground color because\n");
          fprintf(new_fp,"#flip_backgrounds.sh can change the maps' text from\n");
          fprintf(new_fp,"#black to white, and their backgrounds from white to black.\n");
          fprintf(new_fp,"$gmt_prefix gmtset COLOR_BACKGROUND=2/2/2 COLOR_FOREGROUND=253/253/253\n");
          fprintf(new_fp,"digits=%d\n",plot_options.scale_digits);
          fprintf(new_fp,"$gmt_prefix gmtset D_FORMAT=%%.${digits}f\n");
        }
        //Define the filenames for PostScript output.
        sprintf(s,"%s_%04lld",output_base.c_str(),i+1);
        fprintf(new_fp,"#######################################################\n");
        fprintf(new_fp,"data_name=\"%s\"\n",plot_options.just_the_filenames.at(i).c_str());
        fprintf(new_fp,"plot_base=\"%s\"\n",s);
        fprintf(new_fp,"let index=$index+1\n");
        fprintf(new_fp,"#######################################################\n");
        //Record the filename bases in the flip_backgrounds.sh and trim.sh scripts.
        if(i==0){
          if(plot_options.just_the_filenames.size()==1){
            fprintf(flip_fp,"all_bases=\"%s\"\n",s);//First and last base starts and ends the string list.
            fprintf(trim_fp,"all_bases=\"%s\"\n",s);//First and last base starts and ends the string list.
          }
          else{
            fprintf(flip_fp,"all_bases=\"%s\n",s);//First base starts the string list.
            fprintf(trim_fp,"all_bases=\"%s\n",s);//First base starts the string list.
          }
        }
        else if(i<(long long)plot_options.just_the_filenames.size()-1){
          fprintf(flip_fp,"%s\n",s);//Bases that aren't first or last are unadorned.
          fprintf(trim_fp,"%s\n",s);//Bases that aren't first or last are unadorned.
        }
        else{
          fprintf(flip_fp,"%s\"\n",s);//Last base ends with a " to terminate string list.
          fprintf(trim_fp,"%s\"\n",s);//Last base ends with a " to terminate string list.
        }

        //Plot data, with THE SAME title on top and coastlines.
        write_gmt_map_data(results,grid,plot_options,new_fp,results.titles.at(0),i);

        fprintf(new_fp,"ps2raster $png_options $plot_base.ps #Convert PS to PNG format.\n");
        fprintf(new_fp,"#ps2raster -P -Tf $plot_base.ps #Convert PS to PDF, if uncommented.\n");

        //Delete PS file because it's twice as large (more at -E2000 and 0.25x0.25 global- 2.7MB PDF, 165MB PS!) as the PDF or PNG.
        fprintf(new_fp,"rm -f $plot_base.ps\n");
      }
      break;
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    case 5:{//Map using latlon format.

      output_base = "map_parameter";
      //If RGB, only need GMT commands for a single map.
      if(results.rgb.size()==3 and results.latlon.outputs.size()==3) plot_options.just_the_filenames.resize(1);
      for(i=0;i<(long long)plot_options.just_the_filenames.size();i++){
        //If this is the first file...
        if(i==0){
          write_gmt_defs(new_fp);
          fprintf(new_fp,"color_scheme=%d #1/2=white/black background\n",plot_options.color_scheme);
          fprintf(new_fp,"montage=%d #1=left-justify titles, add (a),(b), run montage.sh.\n",plot_options.montage);
          fprintf(new_fp,"prefixes=('(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' '(k)' '(l)' '(m)' '(n)' '(o)' '(p)' '(q)' '(r)' '(s)' '(t)' '(u)' '(v)' '(w)' '(x)' '(y)' '(z)')\n");
          fprintf(new_fp,"index=-1 #Increments on each map, accesses prefixes above for montage.\n");
          fprintf(new_fp,"png_options=\" -P -Tg \" #PDF default: -E720, else 300 dpi.\n");
          fprintf(new_fp,"if [ $montage != 0 ]\n");
          fprintf(new_fp,"then\n");
          fprintf(new_fp,"  png_options=\" -A\"$png_options\n");
          fprintf(new_fp,"fi\n");
          fprintf(new_fp,"#Force off-white(dark gray) fore(back)ground color because\n");
          fprintf(new_fp,"#flip_backgrounds.sh can change the maps' text from\n");
          fprintf(new_fp,"#black to white, and their backgrounds from white to black.\n");
          fprintf(new_fp,"$gmt_prefix gmtset COLOR_BACKGROUND=2/2/2 COLOR_FOREGROUND=253/253/253\n");
          fprintf(new_fp,"digits=%d\n",plot_options.scale_digits);
          fprintf(new_fp,"$gmt_prefix gmtset D_FORMAT=%%.${digits}f\n");
        }
        //Define the filenames for PostScript output.
        sprintf(s,"%s_%04lld",output_base.c_str(),i+1);
        fprintf(new_fp,"#######################################################\n");
        if(results.rgb.size()==3 and results.latlon.outputs.size()==3) fprintf(new_fp,"data_name=redgreenblue\n");
        else fprintf(new_fp,"data_name=\"%s\"\n",plot_options.just_the_filenames.at(i).c_str());
        fprintf(new_fp,"plot_base=\"%s\"\n",s);
        fprintf(new_fp,"let index=$index+1\n");
        fprintf(new_fp,"#######################################################\n");
        //Record the filename bases in the flip_backgrounds.sh and trim.sh scripts.
        if(i==0){
          if(plot_options.just_the_filenames.size()==1){
            fprintf(flip_fp,"all_bases=\"%s\"\n",s);//First and last base starts and ends the string list.
            fprintf(trim_fp,"all_bases=\"%s\"\n",s);//First and last base starts and ends the string list.
          }
          else{
            fprintf(flip_fp,"all_bases=\"%s\n",s);//First base starts the string list.
            fprintf(trim_fp,"all_bases=\"%s\n",s);//First base starts the string list.
          }
        }
        else if(i<(long long)plot_options.just_the_filenames.size()-1){
          fprintf(flip_fp,"%s\n",s);//Bases that aren't first or last are unadorned.
          fprintf(trim_fp,"%s\n",s);//Bases that aren't first or last are unadorned.
        }
        else{
          fprintf(flip_fp,"%s\"\n",s);//Last base ends with a " to terminate string list.
          fprintf(trim_fp,"%s\"\n",s);//Last base ends with a " to terminate string list.
        }

        //Plot data, with title on top and coastlines.
        write_gmt_map_data(results,grid,plot_options,new_fp,results.titles.at(i),i);

        fprintf(new_fp,"ps2raster $png_options $plot_base.ps #Convert PS to PNG format.\n");
        fprintf(new_fp,"#ps2raster -P -Tf $plot_base.ps #Convert PS to PDF, if uncommented.\n");

        //Delete PS file because it's twice as large (more at -E2000 and 0.25x0.25 global- 2.7MB PDF, 165MB PS!) as the PDF or PNG.
        fprintf(new_fp,"rm -f $plot_base.ps\n");
        
        //Move cpt file to backup version so it won't be automatically used
        //if this script is executed again... just in case write_gmt_cpt fails:
        //it's good for that failure to be obvious.
        if(results.rgb.size() != 3) fprintf(new_fp,"mv map.cpt Zbackup_cpt_$plot_base.cpt\n");
        
        //When this plot is finished, store its data_name in previous_data_name
        //so that the next plot (if it's a phase plot with amplitude masking) can
        //access that data.
        fprintf(new_fp,"previous_data_name=$data_name\n");
      }
      break;
    }
    default: cout<<"!!!!WARNING!!!!!! results.options.output_choice "<<results.options.output_choice<<" isn't recognized."<<endl;
  }//End of switch-case.
  
  //Trim images.
  fprintf(new_fp,"#######################################################\n");
  fprintf(new_fp,". ./trim.sh\n");

  //If requested, change colors using flip_backgrounds.sh.
  fprintf(new_fp,"#######################################################\n");
  fprintf(new_fp,"if [ $color_scheme == 2 ]\n");
  fprintf(new_fp,"then\n");
  fprintf(new_fp,"  . ./flip_backgrounds.sh\n");
  fprintf(new_fp,"fi\n");

  //If requested, make a montage using montage.sh.
  fprintf(new_fp,"#######################################################\n");
  fprintf(new_fp,"if [ $montage != 0 ]\n");
  fprintf(new_fp,"then\n");
  fprintf(new_fp,"  . ./montage.sh\n");
  fprintf(new_fp,"fi\n");
  fclose(new_fp);

  //Either way, finish writing flip_backgrounds.sh so it can be used later.
  finish_flip_backgrounds(flip_fp);
  fclose(flip_fp);
  //Finish writing trim.sh.
  finish_trim(trim_fp);
  fclose(trim_fp);
  
  //Create animate script in outputfolder, which should be netcdf_output.
  extra_file = plot_options.outputfolder;
  extra_file.append("animate.sh");
  extra_fp = fopen(extra_file.c_str(),"w");
  //Warn if file doesn't open correctly.
  if(!extra_fp) cout<<"The animate.sh script couldn't be created."<<endl;
  //Every bash script needs this to be the first line.
  fprintf(extra_fp,"#!/bin/bash\n");
  fprintf(extra_fp,"set -x\n");
  fprintf(extra_fp,"delay=100 #delay in hundredths of a second.\n");
  fprintf(extra_fp,"#size=\"640x480\"\n");
  fprintf(extra_fp,"#size=\"800x600\"\n");
  fprintf(extra_fp,"size=\"1024x768\"\n");
  fprintf(extra_fp,"output_base=\"%s\"\n",output_base.c_str());
  fprintf(extra_fp,"#Imagemagick can also output .mng (animated PNG, not well-supported), but ffmpeg is needed as a delegate for .mp4.\n");
  fprintf(extra_fp,"#convert -verbose -delay $delay -loop 0 $output_base* -resize $size animation.gif\n");
  fprintf(extra_fp,"#Or ffmpeg can output .mp4 directly.\n");
  fprintf(extra_fp,"ffmpeg -f image2 -i $output_base%%d.png animation.mp4\n");
  fclose(extra_fp);

  //Create montage script in outputfolder, which should be netcdf_output.
  extra_file = plot_options.outputfolder;
  extra_file.append("montage.sh");
  extra_fp = fopen(extra_file.c_str(),"w");
  //Warn if file doesn't open correctly.
  if(!extra_fp) cout<<"The montage.sh script couldn't be created."<<endl;
  //Every bash script needs this to be the first line.
  fprintf(extra_fp,"#!/bin/bash\n");
  fprintf(extra_fp,"set -x\n");
  fprintf(extra_fp,"output_base=\"%s\"\n",output_base.c_str());
  fprintf(extra_fp,"montage $output_base* -geometry +2+2 montage.png\n");
  fclose(extra_fp);
}

void run_gmt_script(plot_options_s &plot_options){
  /**********************************************************************
  Purpose: This function runs the GMT script create_plots.sh if requested.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  string temp_string;
  FILE *output;

  //Make GMT script executable and then run it (run only if requested).
  time_t start_time = time(NULL);
  temp_string = "cd \"";
  temp_string.append(plot_options.outputfolder);
  temp_string.append("\" ; chmod +x create_plots.sh ; chmod +x flip_backgrounds.sh ; chmod +x animate.sh ; chmod +x montage.sh ; chmod +x trim.sh ; cp create_plots.sh Zbackup_create_plots.sh ; cp projections.sh Zbackup_projections.sh");
  if(plot_options.no_gmt_plots == 0){
    temp_string.append(" ; ./create_plots.sh");
    cout<<"GMT is creating PostScript output..."<<endl;
  }
  else cout<<"Saving GMT script create_plots.sh for later use..."<<endl;
  //cout<<"GMT popen string: "<<temp_string<<endl;
  output = popen(temp_string.c_str(),"r");
  pclose(output);
  cout<<"Finished."<<endl;
  time_t end_time = time(NULL);
  if(end_time-start_time > 0) cout<<"GMT script took"<<seconds_to_string(end_time-start_time)<<endl;
  cout<<"!!!!READ: http://www.imagemagick.org/Usage/compose/#copyopacity"<<endl;
}


int write_netcdf_output(results_s &results,
                        grid_s &grid,
                        plot_options_s &plot_options,
                        int final){
  /**********************************************************************
  Purpose: This function takes structured or 2D output, converts it to
          NetCDF output files in a netcdf_output subfolder that are read
          by the GMT script create_plots.sh to create maps/plots.
          It then creates and runs create_plots.sh.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          int final - Set to 0 if more files will be written, or 1
                      to finalize and write the GMT scripts.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          NetCDF files are created which can be read by the GMT script
          create_plots.sh
  **********************************************************************/
  long long j,l,numfiles=0;//counters.
  char s[max_length];
 
  plot_options_s netcdf_plot_options = plot_options;
  
  //Make a new subfolder called netcdf_output and save NetCDF files in it.
  netcdf_plot_options.outputfolder.append("netcdf_output/");
  create_folder(netcdf_plot_options.outputfolder);
  //Copy options.txt from outputfolder (not from "c" because it might have
  //been overwritten!) to newly-created subfolder.
  copy_file("options.txt",plot_options.outputfolder,netcdf_plot_options.outputfolder);
  //Copy scripts from "c".
  string c_folder(gaiafolder());
  c_folder.append("c/");
  copy_file("projections.sh",c_folder,netcdf_plot_options.outputfolder);
  copy_file("overflow.sh",c_folder,netcdf_plot_options.outputfolder);
  copy_file("notation.sh",c_folder,netcdf_plot_options.outputfolder);
  cout<<"Writing NetCDF output to: "<<netcdf_plot_options.outputfolder<<endl;

  //Requires NetCDF.
  #ifdef NETCDF_HERE
  switch(results.options.output_choice){
    case 5:{//Map using latlon format.
      numfiles = results.latlon.outputs.size();
      init_netcdf_filenames(numfiles, netcdf_plot_options, results);
     
      //Change the error behavior of the netCDF C++ API by creating an
      //NcError object. Until it is destroyed, this NcError object will
      //ensure that the netCDF C++ API returns error codes on any
      //failure, prints an error message, and leaves any other error
      //handling to the calling program. In the case of this example, we
      //just exit with an NC_ERR error code.
      NcError err(NcError::verbose_nonfatal);

      long long NLAT = (long long)results.latlon.lat.size();
      long long NLON = (long long)results.latlon.lon.size();

      if(numfiles>0) for(l=0;l<numfiles;l++){
        //Create the file.
        NcFile dataFile(netcdf_plot_options.output_files.at(l).c_str(), NcFile::Replace);

        //Check to see if the file was created.
        if(!dataFile.is_valid()) return NC_ERR;

        //Define the dimensions. NetCDF will hand back an ncDim object for
        //each.
        NcDim *latDim, *lonDim;
        if(!(latDim = dataFile.add_dim("latitude", NLAT))) return NC_ERR;
        if(!(lonDim = dataFile.add_dim("longitude", NLON))) return NC_ERR;
            
        //Define the coordinate variables.
        NcVar *latVar, *lonVar;
        if(!(latVar = dataFile.add_var("latitude", ncDouble, latDim))) return NC_ERR;
        if(!(lonVar = dataFile.add_var("longitude", ncDouble, lonDim))) return NC_ERR;
            
        //Define units attributes for coordinate vars. This attaches a
        //text attribute to each of the coordinate variables, containing
        //the units.
        if(!latVar->add_att("units", "degrees_north")) return NC_ERR;
        if(!lonVar->add_att("units", "degrees_east")) return NC_ERR;
            
        //Define the netCDF variable for the output and elevations.
        NcVar *outputVar, *elevVar;
        if(!(outputVar = dataFile.add_var("output", ncDouble, latDim, lonDim))) return NC_ERR;
        if(!(elevVar = dataFile.add_var("elev", ncDouble, latDim, lonDim))) return NC_ERR;

        //Store title as a global attribute, and redundantly store in long_name 
        //for output variable because the CF conventions "strongly recommend" its use.
        if(!dataFile.add_att("title", results.titles.at(l).c_str())) return NC_ERR;
        if(!outputVar->add_att("long_name", results.titles.at(l).c_str())) return NC_ERR;
        if(!elevVar->add_att("long_name", "Elevation")) return NC_ERR;

        //Define units attributes for data variables.
        if(!outputVar->add_att("units", results.units.at(l).c_str())) return NC_ERR;
        if(!elevVar->add_att("units", "meters")) return NC_ERR;

        //Define phase format attributes for output.
        if(!outputVar->add_att("phase_format", results.options.current_phase_format)) return NC_ERR;

        //Only bother if error bars are right length, and this one is > 0.0.
        if(results.latlon.outputs.size() == results.error_bars.size() and !results.latlon.outputs.empty()){
          if(results.error_bars.at(l) > 0.0){
            if(!outputVar->add_att("error_bar", results.error_bars.at(l))) return NC_ERR;
          }
        }

        //Only record number of input pts if there's only one region.
        //This means that it's either a universal region or a monthly map, in which case
        //details of the input data used apply to the entire map so they should be recorded.
        if(results.options.input_size.size() == results.titles.size()){
          if(!outputVar->add_att("input_points_used", (long)results.options.input_size.at(l))) return NC_ERR;
        }

        //Only record start_time and end_time if start_times vector has >= 1 entry.
        if((long long)results.options.start_times.size() >= 1){
          if(!outputVar->add_att("start_time_UTC_J2000", (long)results.options.start_times.at(l))) return NC_ERR;
          if(!outputVar->add_att("end_time_UTC_J2000", (long)results.options.end_times.at(l))) return NC_ERR;
        }

        //If marker_lats/lons vectors of vectors are right length, record their values.
        if(results.marker_lats.size() == results.latlon.outputs.size() and results.marker_lons.size() == results.latlon.outputs.size() and !results.latlon.outputs.empty()){
          if(!results.marker_lats.at(l).empty() and !results.marker_lons.at(l).empty()){
            //Define the markers dimension.
            NcDim *markersDim;
            if(!(markersDim = dataFile.add_dim("markers", (long)results.marker_lats.at(l).size()))) return NC_ERR;
            //Define the marker_lats/lons variables.
            NcVar *marker_latsVar, *marker_lonsVar;
            if(!(marker_latsVar = dataFile.add_var("marker_lats", ncDouble, markersDim))) return NC_ERR;
            if(!(marker_lonsVar = dataFile.add_var("marker_lons", ncDouble, markersDim))) return NC_ERR;
            //Write the marker_lats/lons variable data to the file.
            if(!marker_latsVar->put(&results.marker_lats.at(l)[0],(long)results.marker_lats.at(l).size())) return NC_ERR;
            if(!marker_lonsVar->put(&results.marker_lons.at(l)[0],(long)results.marker_lons.at(l).size())) return NC_ERR;
          }
        }

        //If mascon_lats/lons vectors aren't empty, record their values.
        if(!results.latlon.mascon_lats.empty()){
          //Define the mascons dimension.
          NcDim *masconsDim;
          if(!(masconsDim = dataFile.add_dim("mascons", (long)results.latlon.mascon_lats.size()))) return NC_ERR;
          //Define the mascon_lats/lons variables.
          NcVar *mascon_latsVar, *mascon_lonsVar;
          if(!(mascon_latsVar = dataFile.add_var("mascon_lats", ncDouble, masconsDim))) return NC_ERR;
          if(!(mascon_lonsVar = dataFile.add_var("mascon_lons", ncDouble, masconsDim))) return NC_ERR;
          //Write the mascon_lats/lons variable data to the file.
          if(!mascon_latsVar->put(&results.latlon.mascon_lats[0],(long)results.latlon.mascon_lats.size())) return NC_ERR;
          if(!mascon_lonsVar->put(&results.latlon.mascon_lons[0],(long)results.latlon.mascon_lons.size())) return NC_ERR;
        }

        if(!outputVar->add_att("_FillValue", results.latlon.mask)) return NC_ERR;

        //Write the coordinate variable data to the file.
        if(!latVar->put(&results.latlon.lat[0], NLAT)) return NC_ERR;
        if(!lonVar->put(&results.latlon.lon[0], NLON)) return NC_ERR;
                 
        //Put the output values in dataFile, one row of lat at a time.
        for(j=0;j<NLAT;j++){
          if(!outputVar->set_cur(j, 0)) return NC_ERR;
          if(!outputVar->put(&results.latlon.outputs[l][j][0], 1, NLON)) return NC_ERR;       
          if(!elevVar->set_cur(j, 0)) return NC_ERR;
        }
        //The file is automatically closed by the destructor. This frees
        //up any internal netCDF resources associated with the file, and
        //flushes any buffers.
      }
      break;
    }
    default: cout<<"!!!!WARNING!!!!!! results.options.output_choice "<<results.options.output_choice<<" isn't recognized."<<endl;
  }//End of switch-case.
  #else 
    cout<<"write_netcdf_output() requires NetCDF. Search definitions.cpp for NETCDF_HERE."<<endl;
  #endif
  //Now save a file that programs can use to easily read the filenames.
  //Also write the GMT and KMZ scripts, and run the GMT script if requested.
  //(Only do this if this is the "final" call to write_netcdf_output!)
  if(final == 1 and numfiles > 0){
    create_output_filenames(netcdf_plot_options,results);
    write_gmt_scripts(netcdf_plot_options,grid,results);
    //Only make KMZ files for contour maps.
    if(results.options.output_choice == 5) write_kmz_script(netcdf_plot_options,grid,results);
    run_gmt_script(netcdf_plot_options);
  }
  return 0;
}

//====================================================================
//Spectral color pipeline (Hughes and Williams 2010), excerpted from
//the full GAIA functions.cpp.

//====================================================================

void create_synthetic_plot(results_s &results, int choice, long long numpoints,
                           double x0, double delta, double param1, double param2){
  /**********************************************************************
  Purpose: This function creates a 2D plot with synthetic data.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
        int choice      - Controls type of synthetic plot:
                     10 - Sinusoid w period param1 and amp. param2.
                     20 - Gaussian w center param1 and half-width param2.
                     30 - Power law: y = x^param1.
                    100 - Blackbody spectrum at T = param1 Kelvin.
                   1000 - Gaussian noise with std. param1.
        -----------------
        llong numpoints - # x-values if > 0, else uses current x-values.
        double x0,delta - If numpoints > 0, x_values[n]=x0+delta*n.
        double param1/2 - Used to calculate y_values.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  string temp_string;
  char s[max_length];
  vector<double> temp_v_double;//Used to push_back vectors of dbl vectors
  vector<string> temp_str_double;//Used to push_back vectors of str vectors
  vector< vector<double> > temp_vv_double;//push_back vs of vs of dbl vectors!
  int i=0,j=0,k;

  if(numpoints > 0){//Only set up new plots if requested.
    results.xy.titles.push_back("Testing");
    results.xy.x_units.push_back("year");
    results.xy.y_units.push_back("cm");
    results.xy.legends.push_back(temp_str_double);

    //This is a brand new plot, so push_back using the new vv.
    results.xy.x_values.push_back(temp_vv_double);
    results.xy.y_values.push_back(temp_vv_double);
    
    results.xy.x_values.at(i).push_back(temp_v_double);
    results.xy.y_values.at(i).push_back(temp_v_double);

    results.xy.legends.at(i).push_back("Example 1");
    results.xy.x_values.at(i).at(j).resize(numpoints);
    results.xy.y_values.at(i).at(j).resize(numpoints);
    for(k=0;k<numpoints;k++) results.xy.x_values.at(i).at(j).at(k) = x0+delta*k;
  }
  else numpoints = results.xy.x_values.at(i).at(j).size();

  if(choice==10){
    double omega = twoPi/param1;
    for(k=0;k<numpoints;k++) results.xy.y_values.at(i).at(j).at(k) = param2*sin(omega*results.xy.x_values.at(i).at(j).at(k));
  }
  else if(choice==20){
    sprintf(s,"%.0f nm",param1);
    if((long long)results.xy.legends.size() > i) if((long long)results.xy.legends.at(i).size() > j) results.xy.legends.at(i).at(j) = s;
    for(k=0;k<numpoints;k++){
      //Same peak value:
      results.xy.y_values.at(i).at(j).at(k) = exp(-pow(results.xy.x_values.at(i).at(j).at(k)-param1,2)/(2*pow(param2,2)));
      //Normalized:
      //results.xy.y_values.at(i).at(j).at(k) = 1/(param2*sqrtTwoPi)*exp(-pow(results.xy.x_values.at(i).at(j).at(k)-param1,2)/(2*pow(param2,2)));
    }
  }
  else if(choice==30){
    sprintf(s,"Power = %.0f",param1);
    if((long long)results.xy.legends.size() > i) if((long long)results.xy.legends.at(i).size() > j) results.xy.legends.at(i).at(j) = s;
    for(k=0;k<numpoints;k++){
      results.xy.y_values.at(i).at(j).at(k) = pow(results.xy.x_values.at(i).at(j).at(k),param1);
    }
  }
  else if(choice==100){
    //Blackbody colors match this: http://www.vendian.org/mncharity/dir3/blackbody/
    double h = 6.62606957E-34; //Planck constant. Units: J*s.
    double c = 299792458; //Lightspeed. Units: m/s.
    double boltzmann = 1.3806488E-23; //Boltzmann constant. Units: J/K.
    sprintf(s,"%.0f K",param1);
    if((long long)results.xy.legends.size() > i) if((long long)results.xy.legends.at(i).size() > j) results.xy.legends.at(i).at(j) = s;
    for(k=0;k<numpoints;k++){
      double wavelength = results.xy.x_values.at(i).at(j).at(k)*1E-9; //Units: m.
      results.xy.y_values.at(i).at(j).at(k) = (2*h*c*c)/(pow(wavelength,5))/(exp((h*c)/(boltzmann*param1*wavelength))-1);
    }
  }
  else if(choice==1000){
    //Requires GSL libraries.
    #ifdef GSL_HERE
      const gsl_rng_type * R;
      gsl_rng * r;
      gsl_rng_env_setup();
      R = gsl_rng_default;
      r = gsl_rng_alloc(R);
      gsl_rng_set(r,time(0));
      for(k=0;k<numpoints;k++) results.xy.y_values.at(i).at(j).at(k) = gsl_ran_gaussian(r,param1);
    #else 
      cout<<"Random number generation requires GSL. Search definitions.cpp for GSL_HERE."<<endl;
    #endif
  }
  else cout<<"!!!!WARNING!!!! choice "<<choice<<" wasn't recognized in create_synthetic_plot."<<endl;
}

void convert_timeseries_to_rgb(results_s &input, results_s &spectrum,
                               results_s &cie, plot_options_s &plot_options,
                               int choice, int nonmask_pts, int verbose){
  /**********************************************************************
  Purpose: This function converts a timeseries to RGB color.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          int choice - 1 = use L-S periodogram power.
                       2 = interpolate to evenly spaced timeseries,
                           use FFT power.
          int nonmask_pts= 1 if this is the first pt being converted.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  convert_timeseries_to_xyz(input,spectrum,cie,plot_options,choice,nonmask_pts,verbose);
  xyz2rgb(spectrum);
}

void convert_timeseries_to_xyz(results_s &input, results_s &spectrum,
                               results_s &cie, plot_options_s &plot_options,
                               int choice, int nonmask_pts, int verbose){
  /**********************************************************************
  Purpose: This function converts a timeseries to XYZ tristimulus values.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          int choice - 1 = use L-S periodogram power.
                       2 = interpolate to evenly spaced timeseries,
                           use FFT power.
          int nonmask_pts= 1 if this is the first pt being converted.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long i;
  //So that CIE functions are only loaded once.
  if(cie.xy.x_values.empty()) load_cie_functions(cie);
  if(choice == 1) spectrum.options.type=112;//Calculate Lomb-Scargle power of input data.
  else if(choice == 2){
    results_s evenly_spaced = input;
    //Create evenly spaced timeseries.
    double delta = (input.xy.x_values.at(0).at(0).back() - input.xy.x_values.at(0).at(0).at(0))/(double)(input.xy.x_values.at(0).at(0).size()-1);
    for(i=1;i<(long long)evenly_spaced.xy.x_values.at(0).at(0).size();i++){
      evenly_spaced.xy.x_values.at(0).at(0).at(i) = i*delta + evenly_spaced.xy.x_values.at(0).at(0).at(0);
    }
    //Interpolate input data to that evenly spaced timeseries.
    interpolate(input,evenly_spaced);
    spectrum.options.type=103;//Calculate FFT of input data, return power and phase.
  }
  else cout<<"!!!!WARNING!!!! choice "<<choice<<" wasn't recognized in convert_timeseries_to_xyz()."<<endl;
  timeseries_analysis(input, spectrum);
  //If a specified period is <= 0, that frequency limit = min/max.
  minmax_output(spectrum,0,0);//plot index, quiet.
  if(spectrum.options.min_period > 0) spectrum.xy.xmax = 1/spectrum.options.min_period;
  if(spectrum.options.max_period > 0) spectrum.xy.xmin = 1/spectrum.options.max_period;
  interpolate_to_trim_x_axis(spectrum);
  convert_spectrum_from_frequency_to_period(spectrum);
  //Only write colorscale if this is the first nonmask pt being examined.
  if(nonmask_pts==1) write_rgb_colorscale(spectrum,cie,plot_options,1);//1/0=verbose/quiet.
  match_x_axes(spectrum,cie);
  interpolate(spectrum,cie);
  spectrum2xyz(spectrum,cie,verbose);
}

void load_cie_functions(results_s &cie){
  /**********************************************************************
  Purpose: This function loads CIE color matching functions from disk.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Color matching functions from http://cvrl.ioo.ucl.ac.uk/cmfs.htm

  Currently using these: http://cvrl.ioo.ucl.ac.uk/database/text/cmfs/ciexyz31.htm
  References:
  Guild, J. (1931). The colorimetric properties of the spectrum. Philosophical Transactions of the Royal Society of London, A230, 149-187.
  Wright, W. D. (1928). A re-determination of the trichromatic coefficients of the spectral colours. Transactions of the Optical Society, 30, 141-164.

  Consider changing to these: http://cvrl.ioo.ucl.ac.uk/database/text/cmfs/ciexyzjv.htm
  References:
  Judd, D. B. (1951). Report of U.S. Secretariat Committee on Colorimetry and Artificial Daylight. In Proceedings of the Twelfth Session of the CIE, Stockholm (vol. 1, pp. 11). Paris: Bureau Central de la CIE.
  Vos, J. J. (1978). Colorimetric and photometric properties of a 2-deg fundamental observer. Color Research and Application, 3, 125-128.

  Other resources:
  http://www.cie.co.at/main/freepubs.html
  http://www.mathworks.com/matlabcentral/fileexchange/7021-spectral-and-xyz-color-functions/content/colorMatchFcn.m
  http://www.poynton.com/ColorFAQ.html
  http://www.poynton.com/GammaFAQ.html
  http://www.brucelindbloom.com/
  **********************************************************************/
  char s[max_length];
  string file,token;
  FILE *in_fp;
  double temp_double;
  long long temp_long;

  //Empty cie in case it's full already.
  results_s empty;
  cie = empty;

  //Override the usual map output.
  cie.options.output_choice = 2;

  cie.xy.titles.push_back("CIE 1931 color matching functions");
  cie.xy.x_units.push_back("nm");
  cie.xy.y_units.push_back("power");

  //Just 1 plot.
  cie.xy.x_values.resize(1);
  cie.xy.y_values.resize(1);
  cie.xy.legends.resize(1);

  //This plot has 3 lines.
  cie.xy.x_values.at(0).resize(3);
  cie.xy.y_values.at(0).resize(3);
  cie.xy.legends.at(0).resize(3);
  cie.xy.legends.at(0).at(0) = "x";
  cie.xy.legends.at(0).at(1) = "y";
  cie.xy.legends.at(0).at(2) = "z";

  file = gaiafolder();
  file.append("c/sealevel_spectra/");
  //file.append("ciexyzj.txt");//Every 5 nm from 370nm to 770nm.
  //file.append("ciexyz31_1.csv");//Every 1 nm from 360nm to 830nm.
  //file.append("ciexyz31_1_trimmed.csv");//Every 1 nm from 380nm to 760nm.
  file.append("ciexyz31_1_trimmed_420nm_690nm.csv");//Every 1 nm from 420nm to 690nm.
  
  //Open input file for reading.
  in_fp = fopen(file.c_str(),"r");
  //Warn if file doesn't open correctly.
  if(!in_fp) cout<<"The CIE file, "<<file<<", failed to open."<<endl;
 
  while(fgets(s,max_length,in_fp)){
    //Read wavelengths, record for all 3 functions.
    token = strtok(s, ",");
    sscanf(token.c_str(),"%lld", &temp_long);
    cie.xy.x_values.at(0).at(0).push_back((double)temp_long);
    cie.xy.x_values.at(0).at(1).push_back((double)temp_long);
    cie.xy.x_values.at(0).at(2).push_back((double)temp_long);
    //Read x.
    token = strtok(NULL, ",");
    sscanf(token.c_str(),"%lf", &temp_double);
    cie.xy.y_values.at(0).at(0).push_back(temp_double);
    //Read y.
    token = strtok(NULL, ",");
    sscanf(token.c_str(),"%lf", &temp_double);
    cie.xy.y_values.at(0).at(1).push_back(temp_double);
    //Read z.
    token = strtok(NULL, "\n");
    sscanf(token.c_str(),"%lf", &temp_double);
    cie.xy.y_values.at(0).at(2).push_back(temp_double);
  }
  fclose(in_fp);//Close file.
}

void spectrum2xyz(results_s &spectrum, results_s &cie, int verbose){
  /**********************************************************************
  Purpose: This function converts a spectrum to XYZ tristimulus values.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          x,y,z - XYZ tristimulus values.
  **********************************************************************/
  long long i,l;

  //Make sure spectrum has values at the same wavelengths as the xyz
  //color matching functions before calling this function.

  //Zero x,y,z just in case.
  spectrum.xyz.resize(3,0.0);

  //Convolve spectrum with color matching functions to obtain
  //XYZ tristimulus values.
  for(i=0;i<(long long)spectrum.xy.x_values.at(0).at(0).size();i++){
    for(l=0;l<3;l++){
      spectrum.xyz.at(l) += spectrum.xy.y_values.at(0).at(0).at(i)*cie.xy.y_values.at(0).at(l).at(i);
    //cout<<"spectrum: "<<spectrum.xy.y_values.at(0).at(0).at(i)<<endl;
    //cout<<"cie: "<<cie.xy.y_values.at(0).at(l).at(i)<<endl;
    //cout<<"total: "<<spectrum.xyz.at(l)<<endl;
    }
  }
  if(verbose) cout<<"(X,Y,Z) = ("<<scientific<<setprecision(7)<<spectrum.xyz.at(0)<<", "<<spectrum.xyz.at(1)<<", "<<spectrum.xyz.at(2)<<")"<<endl;
}

void xyz2rgb(results_s &spectrum){
  /**********************************************************************
  Purpose: This function converts XYZ tristimulus values to RGB.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          x,y,z - XYZ tristimulus values.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          r,g,b - Linear (NOT gamma corrected) RGB values,
                  where 0,0,0 is black; 1,1,1 is white.
  ***********************************************************************
  Values from Hughes and Williams 2010 and Charles Poynton's ColorFAQ:
  http://www.poynton.com/ColorFAQ.html
  **********************************************************************/
  long long i;

  //Requires GSL libraries.
  #ifdef GSL_HERE
    gsl_vector * xyz = gsl_vector_calloc(3);
    gsl_vector * rgb = gsl_vector_calloc(3);
    gsl_matrix * A = gsl_matrix_alloc(3,3);

    for(i=0;i<3;i++) gsl_vector_set(xyz,i,spectrum.xyz.at(i));

    //A and Ainv are reversed in Hughes and Williams 2010 equations A4 and A5,
    //compared to Charles Poynton's ColorFAQ:
    //http://www.poynton.com/ColorFAQ.html
    //... and compared to Wikipedia:
    //https://en.wikipedia.org/wiki/SRGB
    //So Ainv isn't necessary!
    /*gsl_matrix_set(Ainv, 0, 0, 0.41239081);
    gsl_matrix_set(Ainv, 1, 0, 0.21263903);
    gsl_matrix_set(Ainv, 2, 0, 0.019330821);
    gsl_matrix_set(Ainv, 0, 1, 0.35758433);
    gsl_matrix_set(Ainv, 1, 1, 0.71516866);
    gsl_matrix_set(Ainv, 2, 1, 0.11919473);
    gsl_matrix_set(Ainv, 0, 2, 0.18048081);
    gsl_matrix_set(Ainv, 1, 2, 0.072192319);
    gsl_matrix_set(Ainv, 2, 2, 0.95053222);// */

    gsl_matrix_set(A, 0, 0,  3.2409699);
    gsl_matrix_set(A, 1, 0, -0.96924375);
    gsl_matrix_set(A, 2, 0,  0.055630032);
    gsl_matrix_set(A, 0, 1, -1.5373832);
    gsl_matrix_set(A, 1, 1,  1.8759676);
    gsl_matrix_set(A, 2, 1, -0.20397685);
    gsl_matrix_set(A, 0, 2, -0.49861079);
    gsl_matrix_set(A, 1, 2,  0.041555082);
    gsl_matrix_set(A, 2, 2,  1.0569714);

    //Multiply A*xyz to obtain rgb.
    gsl_blas_dgemv(CblasNoTrans,1.0,A,xyz,0.0,rgb);

    //Extract rgb values from GSL vector so they can be returned.
    spectrum.rgb.resize(3);
    for(i=0;i<3;i++) spectrum.rgb.at(i) = gsl_vector_get(rgb,i);

    //If any values are < 0, add enough white light to make them all positive.
    double min = *min_element(spectrum.rgb.begin(),spectrum.rgb.end());
    if(min<0){
      //Make "min" positive as in paper's appendix, and add 1/255 so the rescale value isn't 0.
      min = -min + 1.0/255.0;
      //This factor rescales the luminance back to its original value.
      double factor = spectrum.xyz.at(1)/(spectrum.xyz.at(1)+min);
      //cout<<"NEGATIVE (r,g,b) = ("<<scientific<<setprecision(7)<<spectrum.rgb.at(0)<<", "<<spectrum.rgb.at(1)<<", "<<spectrum.rgb.at(2)<<")"<<endl;
      for(i=0;i<3;i++) spectrum.rgb.at(i) = (spectrum.rgb.at(i)+min)*factor;
    }

    //Free the space taken up by the various matrices and vectors.
    gsl_matrix_free(A);
    gsl_vector_free(xyz);
    gsl_vector_free(rgb);
  #else 
    cout<<"xyz2rgb requires GSL. Search definitions.cpp for GSL_HERE."<<endl;
  #endif
}

void rescale_rgb(results_s &results){
  /**********************************************************************
  Purpose: This function rescales large RGB values using the largest
          R,G,B value in the whole set of maps.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          r,g,b - Linear RGB values, where 0,0,0 is black; 1,1,1 is white.
  **********************************************************************/
  long long i,j,l;

  double max=0;
  for(i=0;i<3;i++){
    minmax_output(results,i,0);//0=quiet.
    if(results.max > max) max = results.max;
  }
  //Even if max < 1, rescale. Otherwise spectra of very small values
  //(such as GRACE relative accelerations) will always be dark.
  cout<<"Rescaling all RGB values by "<<scientific<<setprecision(10)<<max<<endl;
  if(results.options.output_choice == 1){
    if(results.outputs.size() != 3){
      cout<<"!!!!WARNING!!!!!! rescale_rgb requires output with 3 maps (RGB), not "<<results.outputs.size()<<" maps."<<endl;
      return;
    }
    for(l=0;l<3;l++){
      for(i=0;i<(long long)results.outputs.at(0).size();i++){
        results.outputs.at(l).at(i) /= max;
      }
    }
  }
  else if(results.options.output_choice == 5){
    if(results.latlon.outputs.size() != 3){
      cout<<"!!!!WARNING!!!!!! rescale_rgb requires output with 3 maps (RGB), not "<<results.latlon.outputs.size()<<" maps."<<endl;
      return;
    }
    for(l=0;l<3;l++){
      for(j=0;j<(long long)results.latlon.lat.size();j++){
        for(i=0;i<(long long)results.latlon.lon.size();i++){
          results.latlon.outputs.at(l).at(j).at(i) /= max;
        }
      }
    }
  }
  else cout<<"!!!WARNING!!! rescale_rgb() didn't recognize output_choice "<<results.options.output_choice<<endl;
}

void rescale_single_rgb(results_s &results, int verbose){
  /**********************************************************************
  Purpose: This function rescales large RGB values in a single RGB triplet.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          r,g,b - Linear RGB values, where 0,0,0 is black; 1,1,1 is white.
  **********************************************************************/
  long long i;

  double max = *max_element(results.rgb.begin(),results.rgb.end());
  //Even if max < 1, rescale. Otherwise spectra of very small values
  //(such as GRACE relative accelerations) will always be dark.
  if(verbose) cout<<"Initial (r,g,b)= ("<<scientific<<setprecision(7)<<results.rgb.at(0)<<", "<<results.rgb.at(1)<<", "<<results.rgb.at(2)<<")"<<endl;
  if(verbose) cout<<"Rescaling by "<<scientific<<setprecision(10)<<max<<endl;
  for(i=0;i<3;i++) results.rgb.at(i) /= max;
  if(verbose) cout<<"Linear (r,g,b) = ("<<scientific<<setprecision(7)<<results.rgb.at(0)<<", "<<results.rgb.at(1)<<", "<<results.rgb.at(2)<<")"<<endl;
}

void rescale_single_xyz(results_s &results){
  /**********************************************************************
  Purpose: This function takes an XYZ tristimulus value and scales it so
          it's not larger than 1.0.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          x,y,z - XYZ tristimulus values.
  **********************************************************************/
  long long i;

  double max = *max_element(results.xyz.begin(),results.xyz.end());
  //Even if max < 1, rescale. Otherwise spectra of very small values
  //(such as GRACE relative accelerations) will always be dark.
  cout<<"Rescaling by "<<scientific<<setprecision(10)<<max<<endl;
  for(i=0;i<3;i++) results.xyz.at(i) /= max;
  cout<<"(X,Y,Z) = ("<<scientific<<setprecision(7)<<results.xyz.at(0)<<", "<<results.xyz.at(1)<<", "<<results.xyz.at(2)<<")"<<endl;
}

void gamma_correct(results_s &results){
  /**********************************************************************
  Purpose: This function transforms all RGB values in 3 maps
          to nonlinear RGB, otherwise known as gamma correction.
          This corrects for nonlinear monitor display.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          r,g,b - Nonlinear (gamma corrected) RGB values,
                  where 0,0,0 is black; 1,1,1 is white.
  ***********************************************************************
  Values from Hughes and Williams 2010 and Charles Poynton's ColorFAQ:
  http://www.poynton.com/ColorFAQ.html
  **********************************************************************/
  long long i,j,l;
  results_s spectrum;
  spectrum.rgb.resize(3);

  //Transform to nonlinear RGB, otherwise known as gamma correction.
  if(results.options.output_choice == 1){
    if(results.outputs.size() != 3){
      cout<<"!!!!WARNING!!!!!! gamma_correct requires output with 3 maps (RGB), not "<<results.outputs.size()<<" maps."<<endl;
      return;
    }
    for(i=0;i<(long long)results.outputs.at(0).size();i++){
      for(l=0;l<3;l++) spectrum.rgb.at(l) = results.outputs.at(l).at(i);
      gamma_correct_single_rgb(spectrum,0);//0 = quiet.
      for(l=0;l<3;l++) results.outputs.at(l).at(i) = spectrum.rgb.at(l);
    }
  }
  else if(results.options.output_choice == 5){
    if(results.latlon.outputs.size() != 3){
      cout<<"!!!!WARNING!!!!!! gamma_correct requires output with 3 maps (RGB), not "<<results.latlon.outputs.size()<<" maps."<<endl;
      return;
    }
    for(j=0;j<(long long)results.latlon.lat.size();j++){
      for(i=0;i<(long long)results.latlon.lon.size();i++){
        for(l=0;l<3;l++) spectrum.rgb.at(l) = results.latlon.outputs.at(l).at(j).at(i);
        gamma_correct_single_rgb(spectrum,0);//0 = quiet.
        for(l=0;l<3;l++) results.latlon.outputs.at(l).at(j).at(i) = spectrum.rgb.at(l);
      }
    }
  }
  else cout<<"!!!WARNING!!! gamma_correct() didn't recognize output_choice "<<results.options.output_choice<<endl;
}

void gamma_correct_single_rgb(results_s &results, int verbose){
  /**********************************************************************
  Purpose: This function transforms a single RGB triplet to
          nonlinear RGB, otherwise known as gamma correction.
          This corrects for nonlinear monitor display.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          r,g,b - Nonlinear (gamma corrected) RGB values,
                  where 0,0,0 is black; 1,1,1 is white.
  ***********************************************************************
  Values from Hughes and Williams 2010 and Charles Poynton's ColorFAQ:
  http://www.poynton.com/ColorFAQ.html
  **********************************************************************/
  long long i;

  //Transform to nonlinear RGB, otherwise known as gamma correction.
  double gamma_inv = 0.45;
  double crit = 0.018;//RGB values are gamma corrected differently below and above crit.
  double h =  4.506813168;
  double g = -0.09914989;
  double f =  1.09914989;
  for(i=0;i<3;i++){
    //Typo in Hughes and Williams 2010 equation A7, compared to Charles Poynton's GammaFAQ:
    //http://www.poynton.com/GammaFAQ.html
    if(results.rgb.at(i) <= crit) results.rgb.at(i) *= h;
    else results.rgb.at(i) = f*pow(results.rgb.at(i),gamma_inv) + g;
  }
  if(verbose) cout<<"Nonlinear (r,g,b) = ("<<scientific<<setprecision(7)<<results.rgb.at(0)<<", "<<results.rgb.at(1)<<", "<<results.rgb.at(2)<<")"<<endl;
}

void raise_y_to_power(results_s &results, double power){
  /**********************************************************************
  Purpose: This function takes XYZ tristimulus values and scales all of
          them so that luminance Y is raised to "power" while keeping
          chromaticity (xy in xyY, x = x/(X+Y+Z), same for y) constant.
          This effectively brightens darker spots in maps.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          double power - Luminance Y is effectively raised to this power.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          x,y,z - XYZ tristimulus values.
  **********************************************************************/
  long long i,j,l;
  results_s spectrum;
  spectrum.xyz.resize(3);

  if(results.options.output_choice != 5){
    cout<<"!!!!WARNING!!!!!! raise_y_to_power can't use results.options.output_choice "<<results.options.output_choice<<endl;
    return;
  }
  if(results.latlon.outputs.size() != 3){
    cout<<"!!!!WARNING!!!!!! raise_y_to_power requires latlon output with 3 maps (RGB- actually XYZ here!), not "<<results.latlon.outputs.size()<<" maps."<<endl;
    return;
  }

  for(j=0;j<(long long)results.latlon.lat.size();j++){
    for(i=0;i<(long long)results.latlon.lon.size();i++){
      for(l=0;l<3;l++) spectrum.xyz.at(l) = results.latlon.outputs.at(l).at(j).at(i);
      raise_single_y_to_power(spectrum,power);
      for(l=0;l<3;l++) results.latlon.outputs.at(l).at(j).at(i) = spectrum.xyz.at(l);
    }
  }
}

void raise_single_y_to_power(results_s &results, double power){
  /**********************************************************************
  Purpose: This function takes a SINGLE XYZ tristimulus triplet and scales
          it so that luminance Y is raised to "power" while keeping
          chromaticity (xy in xyY, x = x/(X+Y+Z), same for y) constant.
          This effectively brightens the resulting color.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          double power - Luminance Y is effectively raised to this power.
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          x,y,z - XYZ tristimulus values.
  **********************************************************************/
  long long l;

  //Divide all X,Y,Z by Y^(1-power), which raises Y to "power".
  power = 1 - power;
  double factor = pow(results.xyz.at(1),power);
  for(l=0;l<3;l++) results.xyz.at(l) /= factor;
}

void match_x_axes(results_s &orig, results_s &match){
  /**********************************************************************
  Purpose: Scale "orig" x axis to match "match"; store in output. 
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long i=0,j;
  results_s output = orig;

  output.xy.x_units = match.xy.x_units;//Matching x-axis matches units too.

  double orig_span = orig.xy.x_values.at(i).at(0).back()-orig.xy.x_values.at(i).at(0).at(0);
  double match_span = match.xy.x_values.at(i).at(0).back()-match.xy.x_values.at(i).at(0).at(0);
  for(j=0;j<(long long)orig.xy.x_values.at(i).at(0).size();j++){
    output.xy.x_values.at(i).at(0).at(j) = match.xy.x_values.at(i).at(0).at(0) + (orig.xy.x_values.at(i).at(0).at(j) - orig.xy.x_values.at(i).at(0).at(0))*match_span/orig_span;
  }
  orig = output;
}

void convert_spectrum_from_frequency_to_period(results_s &orig){
  /**********************************************************************
  Purpose: Convert spectrum from frequency (in Hz) to period (in seconds). 
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long i=0,j;
  long long N = orig.xy.x_values.at(i).at(0).size();
  results_s copy = orig;

  if(orig.xy.x_units.at(i) == "Frequency (1/year)") orig.xy.x_units.at(i) = "Period (years)";
  else if(orig.xy.x_units.at(i) == "Frequency (Hz)") orig.xy.x_units.at(i) = "Period (seconds)";
  else cout<<"!!!WARNING!!! convert_spectrum_from_frequency_to_period() didn't recognize units "<<orig.xy.x_units.at(i)<<endl;
  //Need to reverse order of x_values and y_values.
  for(j=0;j<N;j++){
    orig.xy.x_values.at(i).at(0).at(j) = 1/copy.xy.x_values.at(i).at(0).at(N-j-1);
    orig.xy.y_values.at(i).at(0).at(j) = copy.xy.y_values.at(i).at(0).at(N-j-1);
  }
  //Hughes and Williams 2010: Equal power requires I(period) ~ S(frequency)*frequency^2
  //Reversed order of x_values and y_values.
  for(j=0;j<N;j++) orig.xy.y_values.at(i).at(0).at(j) *= pow(copy.xy.x_values.at(i).at(0).at(N-j-1),2);
}

void spectra_init(results_s &spectra, grid_s &grid, int mod_choice, int &choice, int &norm_rgb){
  /**********************************************************************
  Purpose: This function initializes results for spectral maps.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long j;
  char s[max_length];
  input_s dummy_input;

  if(mod_choice == 2003 or mod_choice == 2004 or mod_choice == 2005) norm_rgb = 0;//Over-ride b/c these choices return Y.
  if(mod_choice == 2000 or mod_choice == 2002 or mod_choice == 2003 or mod_choice == 2005) choice = 1;//Use L-S power.
  else if(mod_choice == 2001 or mod_choice == 2004) choice = 2;//Interpolate to evenly spaced timeseries, use FFT power.
  else cout<<"!!!!WARNING!!!! mod_choice "<<mod_choice<<" wasn't recognized."<<endl;
  spectra.units.at(0) = "period (years)";
  spectra.rgb.resize(3);//This tells NetCDF functions to make an RGB map.
  if(spectra.options.output_choice == 1){
    //Now create all the output vectors, including stuff normally done in define_grid.
    create_output_vectors(spectra,dummy_input,1);//1 = GAIA.
    //Now resize the second dimension (have to do this multiple times)
    //so that the second dimension holds all grid points.
    //Note that the first dimension was set in create_output_vectors().
    for(j=0;j<(long long)spectra.options.output_type.size();j++){
      spectra.outputs.at(j).resize(grid.lat.size(),0.0);
    }
  }
  else if(spectra.options.output_choice == 5){
    //Initialize latlon outputs and elevs.
    init_latlon(spectra,3);
  }
  else cout<<"!!!WARNING!!! spectra_init() didn't recognize output_choice "<<spectra.options.output_choice<<endl;
  if(spectra.options.min_period > 0.0 and spectra.options.max_period > 0.0) sprintf(s,"Spectra, %.2f - %.2f year periods",spectra.options.min_period,spectra.options.max_period);
  else sprintf(s,"Spectra, all periods");
  spectra.titles.at(0) = s;//"red" title shows up on RGB map.
  spectra.titles.at(1) = "green";
  spectra.titles.at(2) = "blue";
}

void spectra_loop(results_s &timeseries, results_s &spectra, results_s &cie, results_s &coefs,
                  int &nonmask_pts, int mod_choice, int &choice, int &simple_detrend,
                  int &norm_rgb, double &parameter1, plot_options_s &plot_options){
  /**********************************************************************
  Purpose: This function is called in a loop. Each time it converts a
          timeseries to a single color for a spectral map.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long l;
  results_s spectrum;
  if(mod_choice == 2002 or mod_choice == 2005){
    results_s evenly_spaced = timeseries;
    //Create evenly spaced timeseries.
    double delta = (timeseries.xy.x_values.at(0).at(0).back() - timeseries.xy.x_values.at(0).at(0).at(0))/(double)(timeseries.xy.x_values.at(0).at(0).size()-1);
    for(l=1;l<(long long)evenly_spaced.xy.x_values.at(0).at(0).size();l++){
      evenly_spaced.xy.x_values.at(0).at(0).at(l) = l*delta + evenly_spaced.xy.x_values.at(0).at(0).at(0);
    }
    //Interpolate input data to that evenly spaced timeseries.
    interpolate(timeseries,evenly_spaced);
  }
  if(simple_detrend) detrend_timeseries(timeseries);
  else remove_fit_from_timeseries(timeseries,coefs);
  spectrum.options.min_period = spectra.options.min_period;
  spectrum.options.max_period = spectra.options.max_period;
  if(norm_rgb){
    convert_timeseries_to_rgb(timeseries,spectrum,cie,plot_options,choice,nonmask_pts,0);//0=quiet.
    double current_max = *max_element(spectrum.rgb.begin(),spectrum.rgb.end());
    if(current_max > parameter1) for(l=0;l<3;l++) spectrum.rgb.at(l) /= current_max/parameter1;
    for(l=0;l<3;l++) spectra.rgb.at(l) = spectrum.rgb.at(l);
  }
  else{
    convert_timeseries_to_xyz(timeseries,spectrum,cie,plot_options,choice,nonmask_pts,0);//0=quiet.
    //The 3 maps temporarily hold XYZ values.
    double current_max = *max_element(spectrum.xyz.begin(),spectrum.xyz.end());
    if(current_max > parameter1) for(l=0;l<3;l++) spectrum.xyz.at(l) /= current_max/parameter1;
    for(l=0;l<3;l++) spectra.rgb.at(l) = spectrum.xyz.at(l);
  }
}

void spectra_end(results_s &spectra, int mod_choice, int norm_rgb){
  /**********************************************************************
  Purpose: This function finishes the spectral maps.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long i,j,l;

  if(norm_rgb) rescale_rgb(spectra);
  else{
    //Rescale XYZ values, which are stored just like RGB.
    rescale_rgb(spectra);
    //raise_y_to_power(spectra,0.333);
    if(mod_choice == 2003 or mod_choice == 2004 or mod_choice == 2005){
      spectra.rgb.clear();
      spectra.xyz.clear();
      spectra.titles.at(0) = "X after normalizing";
      spectra.titles.at(1) = "Y (luminance) after normalizing";
      spectra.titles.at(2) = "Z after normalizing";
      for(l=0;l<3;l++) spectra.units.at(l) = "#";
      //Easier to plot if they're 0-255.
      if(spectra.options.output_choice == 1){
        for(i=0;i<(long long)spectra.outputs.at(0).size();i++){
          for(l=0;l<3;l++) spectra.outputs.at(l).at(i) *= 255.0;
        }
      }
      else if(spectra.options.output_choice == 5){
        for(j=0;j<(long long)spectra.latlon.lat.size();j++){
          for(i=0;i<(long long)spectra.latlon.lon.size();i++){
            for(l=0;l<3;l++) spectra.latlon.outputs.at(l).at(j).at(i) *= 255.0;
          }
        }
      }
      else cout<<"!!!WARNING!!! spectra_end() didn't recognize output_choice "<<spectra.options.output_choice<<endl;
      return;
    }
    //Convert XYZ values to RGB.
    results_s spectrum;
    spectrum.xyz.resize(3);
    if(spectra.options.output_choice == 1){
      for(i=0;i<(long long)spectra.outputs.at(0).size();i++){
        for(l=0;l<3;l++) spectrum.xyz.at(l) = spectra.outputs.at(l).at(i);
        xyz2rgb(spectrum);
        for(l=0;l<3;l++) spectra.outputs.at(l).at(i) = spectrum.rgb.at(l);
      }
    }
    else if(spectra.options.output_choice == 5){
      for(j=0;j<(long long)spectra.latlon.lat.size();j++){
        for(i=0;i<(long long)spectra.latlon.lon.size();i++){
          for(l=0;l<3;l++) spectrum.xyz.at(l) = spectra.latlon.outputs.at(l).at(j).at(i);
          xyz2rgb(spectrum);
          for(l=0;l<3;l++) spectra.latlon.outputs.at(l).at(j).at(i) = spectrum.rgb.at(l);
        }
      }
    }
    else cout<<"!!!WARNING!!! spectra_end() didn't recognize output_choice "<<spectra.options.output_choice<<endl;
    rescale_rgb(spectra);//After XYZ is normalized, do the RGB values need normalizing?
  }
  gamma_correct(spectra);
  //Convert RGB values from 0-1 to 0-255 for display.
  if(spectra.options.output_choice == 1){
    for(i=0;i<(long long)spectra.outputs.at(0).size();i++){
      for(l=0;l<3;l++) spectra.outputs.at(l).at(i) *= 255.0;
    }
  }
  else if(spectra.options.output_choice == 5){
    for(j=0;j<(long long)spectra.latlon.lat.size();j++){
      for(i=0;i<(long long)spectra.latlon.lon.size();i++){
        for(l=0;l<3;l++) spectra.latlon.outputs.at(l).at(j).at(i) *= 255.0;
      }
    }
  }
  else cout<<"!!!WARNING!!! spectra_end() didn't recognize output_choice "<<spectra.options.output_choice<<endl;
}

void interpolate(results_s &orig, results_s &match){
  /**********************************************************************
  Purpose: Interpolate x_values of 0th series in 0th plot of orig to 
          x_values of 0th series in 0th plot in match, return as "orig".
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long i=0,j;

  results_s output = orig;

  //Output should have "match" x_values and y_values.
  output.xy.x_values.at(i).at(0) = match.xy.x_values.at(i).at(0);
  output.xy.y_values.at(i).at(0) = match.xy.y_values.at(i).at(0);
  
  //Set y_values to zero so it's obvious if they're not interpolated.
  for(j=0;j<(long long)output.xy.y_values.at(i).at(0).size();j++) output.xy.y_values.at(i).at(0).at(j)=0;

  //Just to make sure interpolation starts with a blank slate.
  orig.xy.closest_indices.clear();
  orig.xy.differences.clear();
  for(j=0;j<(long long)match.xy.x_values.at(i).at(0).size();j++){
    double match_pt = match.xy.x_values.at(i).at(0).at(j);
    if(match_pt >= orig.xy.x_values.at(i).at(0).at(0) and match_pt <= orig.xy.x_values.at(i).at(0).back()){
      closest_points(orig,match, i, j);
      //Use closest_indices to interpolate orig.xy values to match.xy.x_values.
      double x_factor = (match_pt-orig.xy.x_values.at(i).at(0).at(orig.xy.closest_indices.at(1)))/(orig.xy.x_values.at(i).at(0).at(orig.xy.closest_indices.at(1))-orig.xy.x_values.at(i).at(0).at(orig.xy.closest_indices.at(0)));
      //cout<<"x_factor = "<<x_factor<<endl;
      output.xy.y_values.at(i).at(0).at(j) = (orig.xy.y_values.at(i).at(0).at(orig.xy.closest_indices.at(1))-orig.xy.y_values.at(i).at(0).at(orig.xy.closest_indices.at(0)))*x_factor + orig.xy.y_values.at(i).at(0).at(orig.xy.closest_indices.at(1));
    }
    else cout<<"!!!WARNING!!! match_pt "<<scientific<<setprecision(10)<<match_pt<<" outside orig.xy.x_values range: ["<<orig.xy.x_values.at(i).at(0).at(0)<<" , "<<orig.xy.x_values.at(i).at(0).back()<<"]."<<endl;
  }
  //Print first and last nonzero values in this series of y_values.
  /*long long nonzero1=-1,nonzero2=-1;
  for(j=0;j<(long long)output.xy.y_values.at(i).at(0).size();j++){
    if(output.xy.y_values.at(i).at(0).at(j) != 0.0){
      nonzero2 = j;//Record "last" nonzero value regardless.
      if(nonzero1 == -1) nonzero1 = j;//Only if this is the FIRST nonzero value.
    }
  }
  cout<<"After interpolating "<<orig.xy.titles.at(i)<<endl;
  if(nonzero1 != -1) cout<<"First nonzero value at index "<<nonzero1<<" (UTC "<<match.xy.x_values.at(i).at(0).at(nonzero1)<<") is "<<output.xy.y_values.at(i).at(0).at(nonzero1)<<endl;
  else cout<<"No nonzero values at all!"<<endl;
  if(nonzero2 != -1) cout<<"Last nonzero value at index "<<nonzero2<<" (UTC "<<match.xy.x_values.at(i).at(0).at(nonzero2)<<") is "<<output.xy.y_values.at(i).at(0).at(nonzero2)<<endl;// */
  orig = output;
}

void closest_points(results_s &orig, results_s &match, long long i, long long j){
  /**********************************************************************
  Purpose: Find closest points in orig (i) to match point j.
          If closest_indices has size num_indices,
          assume point "j-1" was just used in the same manner.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  long long k;
  double temp_double;
  long long temp_long;
  orig.xy.numindices = 2;//Number of closest matching indices.
  if(orig.xy.closest_indices.size() != orig.xy.differences.size()) cout<<"!!!WARNING!!! orig.xy.closest_indices.size() = "<<orig.xy.closest_indices.size()<<" while orig.xy.differences.size() = "<<orig.xy.differences.size()<<endl;
  
  //Most common case- assume point "j-1" was matched previously.
  if(orig.xy.closest_indices.size() == orig.xy.numindices){
    //Find new differences between closest matches to the last match pt.
    for(k=0;k<(long long)orig.xy.closest_indices.size();k++){
      orig.xy.differences.at(k) = orig.xy.x_values.at(i).at(0).at(orig.xy.closest_indices.at(k)) - match.xy.x_values.at(0).at(0).at(j);
    }
    //Do these matches need to be re-ordered?
    for(k=0;k<(long long)orig.xy.closest_indices.size()-1;k++){
      if(fabs(orig.xy.differences.at(k)) > fabs(orig.xy.differences.at(k+1))){
        temp_double = orig.xy.differences.at(k);
        temp_long = orig.xy.closest_indices.at(k);
        orig.xy.differences.at(k)       = orig.xy.differences.at(k+1);
        orig.xy.closest_indices.at(k)   = orig.xy.closest_indices.at(k+1);
        orig.xy.differences.at(k+1)     = temp_double;
        orig.xy.closest_indices.at(k+1) = temp_long;
      }
    }

    orig.xy.max_diff = *max_element(orig.xy.differences.begin(),orig.xy.differences.end());
    orig.xy.min_diff = *min_element(orig.xy.differences.begin(),orig.xy.differences.end());
    orig.xy.largest_index  = *max_element(orig.xy.closest_indices.begin(),orig.xy.closest_indices.end());
    orig.xy.smallest_index = *min_element(orig.xy.closest_indices.begin(),orig.xy.closest_indices.end());

    //If max difference is non-positive, orig x_values are all less than match_pt,
    //so try next index in list (unless it's the last index).
    if(orig.xy.max_diff <= 0 and orig.xy.largest_index < (long long)orig.xy.x_values.at(i).at(0).size()-1){
      while(is_this_point_closer(orig,match,i,j,orig.xy.largest_index+1) and orig.xy.largest_index < (long long)orig.xy.x_values.at(i).at(0).size()-1){ }
    }
    //Else if min difference is non-negative, orig x_values are all greater than match_pt,
    //so try index before earliest index in list (unless the earliest index in list
    //is the first pt).
    else if(orig.xy.min_diff >= 0 and orig.xy.smallest_index > 0){
      while(is_this_point_closer(orig,match,i,j,orig.xy.smallest_index-1)){ }
    }
  }
  //If closest_indices is cleared, start at the beginning of orig series.
  else if(orig.xy.closest_indices.empty()){
    //Start at the beginning of the orig series.
    for(k=0;k<orig.xy.numindices;k++){
      if(is_this_point_closer(orig,match,i,j,k)){ }
    }
    //If max difference is non-positive, orig x_values are all less than match_pt,
    //so try next index in list (unless it's the last index).
    if(orig.xy.max_diff <= 0 and orig.xy.largest_index < (long long)orig.xy.x_values.at(i).at(0).size()-1){
      while(is_this_point_closer(orig,match,i,j,orig.xy.largest_index+1) and orig.xy.largest_index < (long long)orig.xy.x_values.at(i).at(0).size()-1){ }
    }
  }
  else{
    cout<<"!!!WARNING!!! orig.xy.closest_indices.size() = "<<orig.xy.closest_indices.size()<<endl;
    cout<<"!!!WARNING!!! orig.xy.differences.size() = "<<orig.xy.differences.size()<<endl;
  }
  //cout<<"match pt "<<j<<" x_value: "<<match.xy.x_values.at(0).at(0).at(j)<<endl;
  //for(k=0;k<(long long)orig.xy.closest_indices.size();k++){
    //cout<<"closest indices["<<k<<"] = "<<orig.xy.closest_indices.at(k)<<" and x_value = "<<orig.xy.x_values.at(i).at(0).at(orig.xy.closest_indices.at(k))<<endl;
  //}
}

int is_this_point_closer(results_s &orig, results_s &match,
                         long long i, long long j, long long k){
  /**********************************************************************
  Purpose: Is pt "k" in ith orig series closer to match pt "j"
          than any pts recorded in closest_indices?
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
          Returns "1" if it is closer, "0" if not. Also records 
          closest index and differences.
  **********************************************************************/
  if(orig.xy.closest_indices.size() != orig.xy.differences.size()) cout<<"!!!WARNING!!! orig.xy.closest_indices.size() = "<<orig.xy.closest_indices.size()<<" while orig.xy.differences.size() = "<<orig.xy.differences.size()<<endl;
  if(orig.xy.numindices == 0) cout<<"!!!WARNING!!! orig.xy.numindices = "<<orig.xy.numindices<<endl;

  long long l;
  double diff = orig.xy.x_values.at(i).at(0).at(k) - match.xy.x_values.at(i).at(0).at(j);
  int closer = 0;//Assume this pt ISN'T closer until proven wrong.
  
  //Most common case- already found numindices closest matches.
  if(orig.xy.closest_indices.size() == orig.xy.numindices){
    l=0;//Compare to closest match first, then next-closest, etc.
    while(l<(long long)orig.xy.closest_indices.size() and closer == 0){
      if(fabs(diff) < fabs(orig.xy.differences.at(l))){
        orig.xy.closest_indices.insert(orig.xy.closest_indices.begin()+l,k);
        orig.xy.differences.insert(orig.xy.differences.begin()+l,diff);
        //Now delete farthest match to keep orig.xy.numindices matches only.
        orig.xy.closest_indices.pop_back();
        orig.xy.differences.pop_back();
        //Exit while loop.
        closer = 1;
      }
      else l++;//Compare to next-closest match.
    }
  }
  else if(orig.xy.closest_indices.empty()){//Start from scratch.
    orig.xy.closest_indices.push_back(k);
    orig.xy.differences.push_back(diff);
    closer = 1;//When closest_indices is cleared, this is true by default.
  }
  else if(orig.xy.closest_indices.size() < orig.xy.numindices){
    l=0;//Compare to closest match first, then next-closest, etc.
    while(l<(long long)orig.xy.closest_indices.size() and closer == 0){
      if(fabs(diff) < fabs(orig.xy.differences.at(l))){
        orig.xy.closest_indices.insert(orig.xy.closest_indices.begin()+l,k);
        orig.xy.differences.insert(orig.xy.differences.begin()+l,diff);
        //Exit while loop.
        closer = 1;
      }
      else l++;//Compare to next-closest match.
    }
    if(closer == 0){//Match isn't closer, but closest_indices isn't full.
      orig.xy.closest_indices.push_back(k);
      orig.xy.differences.push_back(diff);
      closer = 1;//Match found- albeit only to fill closest_indices.
    }
  }
  else{
    cout<<"!!!WARNING!!! orig.xy.closest_indices.size() = "<<orig.xy.closest_indices.size()<<endl;
    cout<<"!!!WARNING!!! orig.xy.differences.size() = "<<orig.xy.differences.size()<<endl;
  }
  orig.xy.max_diff = *max_element(orig.xy.differences.begin(),orig.xy.differences.end());
  orig.xy.min_diff = *min_element(orig.xy.differences.begin(),orig.xy.differences.end());
  orig.xy.largest_index  = *max_element(orig.xy.closest_indices.begin(),orig.xy.closest_indices.end());
  orig.xy.smallest_index = *min_element(orig.xy.closest_indices.begin(),orig.xy.closest_indices.end());
  return closer;
}

void finish_flip_backgrounds(FILE *flip_fp){
  /**********************************************************************
  Purpose: This function finishes writing the flip_backgrounds.sh script.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          flip_fp - file pointer to flip_backgrounds.sh script file
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  fprintf(flip_fp,"for current_base in $all_bases\n");
  fprintf(flip_fp,"do\n");
  fprintf(flip_fp,"  #Change black to cyan temporarily.\n");
  fprintf(flip_fp,"  convert $current_base.png -fill cyan -opaque black $current_base.png\n");
  fprintf(flip_fp,"  #Change white to black.\n");
  fprintf(flip_fp,"  convert $current_base.png -fill black -opaque white $current_base.png\n");
  fprintf(flip_fp,"  #Change temporary cyan to white.\n");
  fprintf(flip_fp,"  convert $current_base.png -fill white -opaque cyan $current_base.png\n");
  fprintf(flip_fp,"done\n");
}

void finish_trim(FILE *trim_fp){
  /**********************************************************************
  Purpose: This function finishes writing the trim.sh script.
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
          trim_fp - file pointer to trim.sh script file
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)
  **********************************************************************/
  fprintf(trim_fp,"for current_base in $all_bases\n");
  fprintf(trim_fp,"do\n");
  fprintf(trim_fp,"  convert $current_base.png -trim $current_base.png\n");
  fprintf(trim_fp,"done\n");
}

string gaiafolder(){
  /**********************************************************************
  Purpose: This function assumes that it's being run from
          GAIAFOLDER/c, so it reports the "realpath" of ".."
  ***********************************************************************
  Input:  (types ending in "_s" are defined in definitions.hpp)
  ***********************************************************************
  Output: (types ending in "_s" are defined in definitions.hpp)

  Future upgrades should consider:
  http://insanecoding.blogspot.com/2007/11/implementing-realpath-in-c.html
  http://insanecoding.blogspot.com/2007/11/pathmax-simply-isnt.html
  **********************************************************************/
  char resolved_path[PATH_MAX];
  if(realpath("..", resolved_path)==NULL) cout<<"!!!WARNING!!! realpath() error!"<<endl;
  string output = resolved_path;
  output.append("/");//Folders always end in "/".
  return output;
}
