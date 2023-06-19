import os
import numpy as np
import datetime

from typing import Dict

plot_options = {
    'outputfolder': '.',
    'just_the_filenames': ['r.nc','g.nc','b.nc'],
    'output_base': 'map_parameter',
    'GMT5': 1,
    'projection': 1, # Since the function uses plot_options.projection, a reasonable default might be 1 for the Robinson projection.
    'plot_mascons': 0,
    'coastlines': 1, #1:pscoast, 2:pscoast+InSAR.
    'blurb_disabled': 1, 
    'montage': 0, #1=left-justify titles, add (a),(b), run montage.sh.
    'region': 'global',    # Define a global region for all plots
    'frame': 'a',    # Add a frame to the plot with automatic tick intervals
    'color_scheme': 2, #1/2=white/black background
    'land_color': 'white',   # Fill continents with color white
    'sea_color': 'blue',   # Fill oceans with color blue
    'scale_digits': 2,    # Set the number of digits for the D_FORMAT
    'show_fig': False,
    'save_fig': True,
    'figsize': (10, 10),
    'dpi': 300,
    'linewidth': 2.0,
    'linestyle': '-',
    'color': 'black',
    'marker': 'o',
    'markersize': 5,
    'markerfacecolor': 'blue',
    'markeredgewidth': 1.5,
    'markeredgecolor': 'black',
    'font_size': 14,
    'font_weight': 'bold',
    'x_label': 'X-axis',
    'y_label': 'Y-axis',
    'title': 'My Plot',
    'grid': True,
    'grid_linestyle': '--',
    'grid_linewidth': 0.5,
    'grid_alpha': 0.7
}

grid = {
    'add_on': {
        'write_gmt_defs': None,
        'just_the_filenames': [],
        'color_scheme': 1,  # 1/2=white/black background
        'montage': 0,  # 1=left-justify titles, add (a),(b), run montage.sh.
        'prefixes': ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)', '(q)', '(r)', '(s)', '(t)', '(u)', '(v)', '(w)', '(x)', '(y)', '(z)'],
        'index': -1,  # Increments on each map, accesses prefixes above for montage.
        'png_options': ' -P -Tg ',  # PDF default: -E720, else 300 dpi.
        'digits': 2,  # assumed reasonable default for scale_digits.
        'output_base': '',  # empty string as placeholder, will be updated in each iteration.
        'data_name': '',  # empty string as placeholder, will be updated in each iteration.
        'write_gmt_map_data': None,
        'ps2raster': None,
        'finish_flip_backgrounds': None,
        'finish_trim': None,
        'outputfolder': '',  # empty string as placeholder, need to be updated with real path
        'animate.sh': None,
        'montage.sh': None,
        'write_gmt_colorscale': None,
        'write_rgb_colorscale': None
    }
}

results = {
    'max_widths':2,
    'options': {'output_choice': 5},
    'latlon': {
        'lat': [10.0, 20.0, 30.0, 40.0, 50.0],
        'lon': [10.0, 20.0, 30.0, 40.0, 50.0],
        'outputs': ['output1', 'output2', 'output3'],
        'mascon_lats': [10.0, 20.0, 30.0]
    },
    'marker_lats': ['output1', 'output2', 'output3'],
    'marker_lons': ['output1', 'output2', 'output3'],
    'rgb': [1, 2, 3],
    'rgb_choice': 2,
    'units': ['unit1', 'unit2', 'unit3'],
    'maxlat': 90.0,
    'minlat': -90.0,
    'maxlon': 180.0,
    'minlon': -180.0,
    'error_bars': [0.1, 0.2, 0.3],

    # Assuming the `just_the_filenames` variable is related to the data we are dealing with, 
    # we will set it as a list with a single default value as an example. The actual values 
    # should be replaced according to your specific use case.
    'just_the_filenames': ['example_filename'],

    # Based on the code, `color_scheme` and `montage` seems to be some configuration parameters.
    # I'll assume they are integers and default them to 1.
    'color_scheme': 1,
    'montage': 1,

    # `scale_digits` appears to be controlling the formatting of some GMT parameter. 
    # I'll assume it's an integer and default it to 2 (for 2 decimal places).
    'scale_digits': 2,

    # `output_base` is being used to format filenames, I'll default it to a string 'output'.
    'output_base': 'output',

    # `titles` appears to be a list related to the `just_the_filenames`, hence it should have 
    # the same length. We will initialize it with a single default value as an example.
    'titles': ['example_title'],

    # The actual values of the keys `grid`, `plot_options`, and `outputfolder` are not clear from the provided code,
    # however, their existence can be inferred. For `grid` and `plot_options`, I'll use simple string placeholders,
    # and for `outputfolder`, I'll use a generic path.
    'grid': grid,
    'plot_options': plot_options,
    'outputfolder': '/path/to/outputfolder',
    'date': datetime.datetime.now(),
    'config': {
        'x': 0.0,
        'y': 0.0,
        'width': 0.0,
        'height': 0.0,
    },
    'xy': {
        'x_values': [[0.0]],
        'y_values': [[0.0]],
    },
    'xyz': {
        'x_values': [],
        'y_values': [],
        'z_values': [],
    },
    'min_max': {
        'x': [0.0, 0.0],
        'y': [0.0, 0.0],
        'z': [0.0, 0.0],
    },
    'misc': {
        'scale_format': '',
        'overflow': '',
        'scale_pos': '',
        'units_format': '',
        'units_pos': '',
        'misc_range': '',
        'start': '',
        'middle': '',
        'end': '',
        'scale_width': 0.0,
        'scale_x': 0.0,
        'scale_y': 0.0,
        'scale_length': 0.0,
        'plot_base': '',
        'gmt_prefix': '',
    },
}

def write_gmt_scripts(plot_options, grid, results):

    # Create GMT script in outputfolder, which should be netcdf_output.
    new_file = os.path.join(plot_options['outputfolder'], 'create_plots.sh')
    try:
        new_fp = open(new_file, 'w')
    except IOError:
        print("The create_plots.sh GMT script couldn't be created.")

    # Every bash script needs this to be the first line.
    new_fp.write("#!/bin/bash\n")
    new_fp.write("#set -x #Uncomment to echo these commands.\n")

    # Create flip_backgrounds script in outputfolder, which should be netcdf_output.
    flip_file = os.path.join(plot_options['outputfolder'], 'flip_backgrounds.sh')
    try:
        flip_fp = open(flip_file, 'w')
    except IOError:
        print("The flip_backgrounds.sh script couldn't be created.")

    # Every bash script needs this to be the first line.
    flip_fp.write("#!/bin/bash\n")
    flip_fp.write("set -x\n")

    # Create trim script in outputfolder, which should be netcdf_output.
    trim_file = os.path.join(plot_options['outputfolder'], 'trim.sh')
    try:
        trim_fp = open(trim_file, 'w')
    except IOError:
        print("The trim.sh script couldn't be created.")

    # Every bash script needs this to be the first line.
    trim_fp.write("#!/bin/bash\n")
    trim_fp.write("set -x\n")

    # If RGB, only need GMT commands for a single map.
    if 1:#len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
        just_the_filenames = plot_options['just_the_filenames'][:1]

    for i in range(len(just_the_filenames)):
        # If this is the first file...
        if i == 0:
            write_gmt_defs(new_fp)
            new_fp.write(f"color_scheme={plot_options['color_scheme']} #1/2=white/black background\n")
            new_fp.write(f"montage={plot_options['montage']} #1=left-justify titles, add (a),(b), run montage.sh.\n")
            new_fp.write("prefixes=('(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' '(k)' '(l)' '(m)' '(n)' '(o)' '(p)' '(q)' '(r)' '(s)' '(t)' '(u)' '(v)' '(w)' '(x)' '(y)' '(z)')\n")
            new_fp.write("index=-1 #Increments on each map, accesses prefixes above for montage.\n")
            new_fp.write('png_options=" -P -Tg " #PDF default: -E720, else 300 dpi.\n')
            new_fp.write("if [ $montage != 0 ]\nthen\n  png_options=\" -A\"$png_options\nfi\n")
            new_fp.write("#Force off-white(dark gray) fore(back)ground color because\n#flip_backgrounds.sh can change the maps' text from\n#black to white, and their backgrounds from white to black.\n")
            new_fp.write(f"$gmt_prefix gmtset COLOR_BACKGROUND=2/2/2 COLOR_FOREGROUND=253/253/253\n")
            new_fp.write(f"digits={plot_options['scale_digits']}\n")
            new_fp.write(f"$gmt_prefix gmtset D_FORMAT=%.${{digits}}f\n")

        # Define the filenames for PostScript output.
        s = f"{plot_options['output_base']}_{i+1:04d}"
        new_fp.write("#######################################################\n")
        if 1:#len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
            new_fp.write("data_name=redgreenblue\n")
        else:
            new_fp.write(f'data_name="{just_the_filenames[i]}"\n')
        new_fp.write(f'plot_base="{s}"\n')
        new_fp.write("let index=$index+1\n")
        new_fp.write("#######################################################\n")

        # Record the filename bases in the flip_backgrounds.sh and trim.sh scripts.
        if i == 0:
            if len(just_the_filenames) == 1:
                flip_fp.write(f'all_bases="{s}"\n')  # First and last base starts and ends the string list.
                trim_fp.write(f'all_bases="{s}"\n')  # First and last base starts and ends the string list.
            else:
                flip_fp.write(f'all_bases="{s}\n')  # First base starts the string list.
                trim_fp.write(f'all_bases="{s}\n')  # First base starts the string list.
        elif i < len(just_the_filenames) - 1:
            flip_fp.write(f'{s}\n')  # Bases that aren't first or last are unadorned.
            trim_fp.write(f'{s}\n')  # Bases that aren't first or last are unadorned.
        else:
            flip_fp.write(f'{s}"\n')  # Last base ends with a " to terminate string list.
            trim_fp.write(f'{s}"\n')  # Last base ends with a " to terminate string list.

        # Plot data, with title on top and coastlines.
        write_gmt_map_data(results, grid, plot_options, new_fp, results['titles'][i], i)

        new_fp.write(f"ps2raster $png_options $plot_base.ps #Convert PS to PNG format.\n")
        new_fp.write("#ps2raster -P -Tf $plot_base.ps #Convert PS to PDF, if uncommented.\n")

        # Delete PS file because it's twice as large (more at -E2000 and 0.25x0.25 global- 2.7MB PDF, 165MB PS!) as the PDF or PNG.
        new_fp.write("rm -f $plot_base.ps\n")
        
        # Move cpt file to backup version so it won't be automatically used
        # if this script is executed again... just in case write_gmt_cpt fails:
        # it's good for that failure to be obvious.
        if len(results['rgb']) != 3:
            new_fp.write(f"mv map.cpt Zbackup_cpt_$plot_base.cpt\n")

        # When this plot is finished, store its data_name in previous_data_name
        # so that the next plot (if it's a phase plot with amplitude masking) can
        # access that data.
        new_fp.write("previous_data_name=$data_name\n")
    # Trim images.
    new_fp.write("#######################################################\n")
    new_fp.write(". ./trim.sh\n")

    # If requested, change colors using flip_backgrounds.sh.
    new_fp.write("#######################################################\n")
    new_fp.write("if [ $color_scheme == 2 ]\n")
    new_fp.write("then\n")
    new_fp.write("  . ./flip_backgrounds.sh\n")
    new_fp.write("fi\n")

    # If requested, make a montage using montage.sh.
    new_fp.write("#######################################################\n")
    new_fp.write("if [ $montage != 0 ]\n")
    new_fp.write("then\n")
    new_fp.write("  . ./montage.sh\n")
    new_fp.write("fi\n")
    new_fp.close()

    # Either way, finish writing flip_backgrounds.sh so it can be used later.
    finish_flip_backgrounds(flip_fp)
    flip_fp.close()

    # Finish writing trim.sh.
    finish_trim(trim_fp)
    trim_fp.close()

    # Create animate script in outputfolder, which should be netcdf_output.
    extra_file = os.path.join(plot_options['outputfolder'], "animate.sh")
    try:
        extra_fp = open(extra_file, "w")
    except IOError:
        print("The animate.sh script couldn't be created.")

    # Every bash script needs this to be the first line.
    extra_fp.write("#!/bin/bash\n")
    extra_fp.write("set -x\n")
    extra_fp.write("delay=100 #delay in hundredths of a second.\n")
    extra_fp.write("#size=\"640x480\"\n")
    extra_fp.write("#size=\"800x600\"\n")
    extra_fp.write("size=\"1024x768\"\n")
    extra_fp.write(f"output_base=\"{plot_options['output_base']}\"\n")
    extra_fp.write("#Imagemagick can also output .mng (animated PNG, not well-supported), but ffmpeg is needed as a delegate for .mp4.\n")
    extra_fp.write("#convert -verbose -delay $delay -loop 0 $output_base* -resize $size animation.gif\n")
    extra_fp.write("#Or ffmpeg can output .mp4 directly.\n")
    extra_fp.write("ffmpeg -f image2 -i $output_base%d.png animation.mp4\n")
    extra_fp.close()

    # Create montage script in outputfolder, which should be netcdf_output.
    extra_file = os.path.join(plot_options['outputfolder'], "montage.sh")
    try:
        extra_fp = open(extra_file, "w")
    except IOError:
        print("The montage.sh script couldn't be created.")

    # Every bash script needs this to be the first line.
    extra_fp.write("#!/bin/bash\n")
    extra_fp.write("set -x\n")
    extra_fp.write(f"output_base=\"{plot_options['output_base']}\"\n")
    extra_fp.write("montage $output_base* -geometry +2+2 montage.png\n")
    extra_fp.close()

def write_gmt_coastlines(new_fp):
    """
    Purpose: This function writes the GMT coastlines command(s).
    Input:  
        new_fp - file pointer to current gmt script file
    """
    new_fp.write("$gmt_prefix pscoast -W$pscoast_thk/$coast_color $pscoast_res $range $projection $map_pos $middle >> $plot_base.ps\n")
    new_fp.write("if [ $coastlines == 2 ]\n")
    new_fp.write("then\n")
    new_fp.write("  $gmt_prefix psxy -N $coast_file -: -Sc$coast_thk -W$coast_thk/$coast_color $range $projection $map_pos $middle >> $plot_base.ps\n")
    new_fp.write("fi\n")

def is_polar(results: Dict, grid: Dict) -> int:
    """
    Purpose: This function examines results,grid, and returns 0 if global,
            or 1 (2) if it's roughly centered on the north (south) pole.
    Input: results, grid (types are dictionary)
    Output: Return value described above.
    """
    polar = 0  # 0 - global, 1 - north pole, 2 - south pole.
    
    if results['options']['output_choice'] == 1 or results['options']['output_choice'] == 4:
        results['minlat'] = min(grid['lat'])
        results['maxlat'] = max(grid['lat'])
        results['minlon'] = min(grid['lon'])
        results['maxlon'] = max(grid['lon'])
        
    elif results['options']['output_choice'] == 5:
        results['minlat'] = min(results['latlon']['lat'])
        results['maxlat'] = max(results['latlon']['lat'])
        results['minlon'] = min(results['latlon']['lon'])
        results['maxlon'] = max(results['latlon']['lon'])
        
    else:
        print(f"!!!!WARNING!!!!!! results['options']['output_choice'] {results['options']['output_choice']} isn't recognized.")
        
    if results['minlat'] < 0 and results['maxlat'] > 0:
        polar = 0
    elif results['minlat'] > 0 and results['maxlat'] > 0:
        polar = 1
    elif results['minlat'] < 0 and results['maxlat'] < 0:
        polar = 2
    
    return polar


def write_gmt_colorscale(new_fp, results, kml_output):
    """
    Purpose: This function writes the GMT color scale command.
    Input:  
        new_fp - file pointer to current gmt script file
        kml_output - 0 for GMT scale next to plot, 1 for KMZ by itself.
    """
    rgb = {}  # Just to init choice and maxes.

    if len(results['rgb']) != 3:
        new_fp.write("cpt_name=\"-Cmap.cpt \"\n")
    else:
        new_fp.write("cpt_name=\"-C../../rgb00001.cpt \"\n")
    
    if len(results['rgb']) != 3 or results['rgb_choice'] < 2:
        # If this is intended for KMZ output, don't overlay the scale.
        if kml_output:
            # KMZ version writes to a different file.
            new_fp.write("$gmt_prefix psscale $cpt_name -L $scale_format $overflow $scale_pos -A $start > scale_$plot_base.ps\n")
            new_fp.write("# Print units manually, otherwise they're too close to numbers on scale.\n")
            new_fp.write("echo $units_format $scale_units | $gmt_prefix pstext -N $units_pos $misc_range $end >> scale_$plot_base.ps\n")
        else:
            new_fp.write("$gmt_prefix psscale $cpt_name -L $scale_format $overflow $scale_pos -A $middle >> $plot_base.ps\n")
            new_fp.write("# Print units manually, otherwise they're too close to numbers on scale.\n")
            new_fp.write("echo $units_format $scale_units | $gmt_prefix pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n")
    else:
        new_fp.write("numwidths={}\n".format(results['max_widths']))
        new_fp.write("$gmt_prefix gmtset TICK_LENGTH 0.3c\n")
        new_fp.write("scale_width=$(bc <<< \"scale=5; $scale_width / $numwidths\")\n")
        new_fp.write("for (( j=1; j <= $numwidths; j++ ))\n")
        new_fp.write("do\n")
        new_fp.write("  j_string=$(printf '%%05d' $j)\n")
        new_fp.write("  cpt_name=\"-C../rgb$j_string.cpt\"\n")
        new_fp.write("  $gmt_prefix psscale $cpt_name -L $scale_format $overflow $scale_pos -S -A $middle >> $plot_base.ps\n")
        new_fp.write("  # Print units manually, otherwise they're too close to numbers on scale.\n")
        new_fp.write("  echo $units_format $scale_units | $gmt_prefix pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n")
        new_fp.write("  # Move next scale to the right and get rid of the tick marks.\n")
        new_fp.write("  scale_x=$(bc <<< \"scale=5; $scale_x+$scale_width\")\n")
        new_fp.write("  scale_pos=\" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c \"\n")
        new_fp.write("  $gmt_prefix gmtset TICK_LENGTH 0.0\n")
        new_fp.write("done\n")






def write_rgb_colorscale(results, cie, plot_options, verbose):
    # Initialize objects
    i = j = l = 0
    s = [None]*max_length
    new_file = ''
    old = copy = all = results_s()

    if results['rgb_choice'] == 0:
        new_file = plot_options['outputfolder'] + "rgb00001.cpt"
        with open(new_file, 'w') as new_fp:
            new_fp.write("# COLOR_MODEL = RGB\n")
            print("Writing RGB colorscale to disk.")
            for i in range(len(results['xy']['x_values'][0][0])):
                copy = results
                create_synthetic_plot(copy,20,0,0,0,results['xy']['x_values'][0][0][i],0.2)
                match_x_axes(copy, cie)
                interpolate(copy, cie)
                spectrum2xyz(copy, cie, verbose)
                xyz2rgb(copy)
                rescale_single_rgb(copy, verbose)
                gamma_correct_single_rgb(copy, verbose)
                for l in range(3): copy.rgb[l] *= 255.0
                if i > 0:
                    new_fp.write(f"{results['xy']['x_values'][0][0][i-1]:12.6e} {int(old.rgb[0]):3d} {int(old.rgb[1]):3d} {int(old.rgb[2]):3d} {results['xy']['x_values'][0][0][i]:12.6e} {int(copy.rgb[0]):3d} {int(copy.rgb[1]):3d} {int(copy.rgb[2]):3d}\n")
                old = copy
    elif results['rgb_choice'] == 1:
        # Create CPT file in outputfolder, which is NOT netcdf_output bc that doesn't exist yet.
        new_file = os.path.join(plot_options['outputfolder'], 'rgb00001.cpt')
        try:
            with open(new_file, 'w') as new_fp:
                new_fp.write("# COLOR_MODEL = RGB\n")
                print("Writing RGB colorscale to disk.")
                
                hw = 2 * (results['xy']['x_values'][0][0][1] - results['xy']['x_values'][0][0][0])
                all['options']['output_choice'] = 5
                all['latlon']['lat'] = [0]*len(results['xy']['x_values'][0][0])
                all['latlon']['lon'] = [0]
                j = 0  # Because there's only 1 lon.
                init_latlon(all,3)
                
                # Loop through periods for colorscale divisions.
                for i in range(len(results['xy']['x_values'][0][0])):
                    # Calculate color for each period.
                    copy = deepcopy(results)
                    create_synthetic_plot(copy,20,0,0,0,results['xy']['x_values'][0][0][i],hw) # choice,#pts,x0,delta,param1,2.
                    match_x_axes(copy,cie)
                    interpolate(copy,cie)
                    spectrum2xyz(copy,cie,verbose)
                    xyz2rgb(copy)
                    for l in range(3):
                        all['latlon']['outputs'][l][i][j] = copy.rgb[l]
                
                spectra_end(all,2000,1) # mod_choice,norm_rgb
                
                for i in range(1,len(results['xy']['x_values'][0][0])):
                    new_fp.write("%12.6e %3d %3d %3d %12.6e %3d %3d %3d\n" % (
                        results['xy']['x_values'][0][0][i-1], int(all['latlon']['outputs'][0][i-1][j]), int(all['latlon']['outputs'][1][i-1][j]), int(all['latlon']['outputs'][2][i-1][j]),
                        results['xy']['x_values'][0][0][i], int(all['latlon']['outputs'][0][i][j]), int(all['latlon']['outputs'][1][i][j]), int(all['latlon']['outputs'][2][i][j])
                    ))
        except Exception as e:
            print(f"The rgb00001.cpt file couldn't be created due to the following error: {e}")

    else: print(f"!!!WARNING!!! Didn't recognize {results['rgb_choice'] = }")

def write_gmt_map_data(results, grid, plot_options, new_fp, title, i):
    """
    Purpose: This function writes the GMT data plotting command.
            If the title is blank, this is assumed to be a KMZ map so
            the map fills the page in the cylindrical equidistant
            projection that works in Google Earth.
    Input:  new_fp - file handle to current gmt script file
            title - string to be plotting at top of map data. If blank,
                    title is suppressed and Google Earth KMZ is output.
            i - (output) index of parameter being mapped.
    """

    # Is this a polar or global plot?
    polar = is_polar(results, grid)  # 0 - global, 1 - north pole, 2 - south pole.

    if i == 0:
        new_fp.write("#Set resolution, coast_file, coast_thickness, and coastlines\n")
        new_fp.write("#on first map only because they should be universal.\n")
        coast_file = plot_options['outputfolder'] + "data/ancillary/Rignot/InSAR_GL_Antarctica.txt"
        new_fp.write(f'coast_file="{coast_file}"\n')

        # Only latlon data determines resolution using delta_lat.
        if results['options']['output_choice'] == 5:
            delta_lat = abs(results['latlon']['lat'][1] - results['latlon']['lat'][0])
            if delta_lat < 0.4:  # Latlon lat spacing controls the resolution.
                new_fp.write('resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write('pscoast_res=" -Df+ "\n')
            elif delta_lat < 0.9:  # Latlon lat spacing controls the resolution.
                new_fp.write('resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write('pscoast_res=" -Df+ "\n')
            else:
                new_fp.write('resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write('pscoast_res=" -Di+ "\n')
        elif results['options']['output_choice'] in [1, 4]:
            new_fp.write('resolution=" -E50 " #50/2000 is low/high quality.\n')
            new_fp.write('pscoast_res=" -Di+ "\n')
        else:
            print(f"!!!!WARNING!!!!!! results['options']['output_choice'] {results['options']['output_choice']} isn't recognized.")

        new_fp.write('pscoast_res_orig=$pscoast_res #Don\'t want USA maps to repeatedly add -N2.\n')
        new_fp.write('pscoast_thk="0.6"\n')
        new_fp.write('coast_thk="0.009"\n')
        if polar == 1:
            coastlines = 1  # InSAR is only in Antarctica, so disable for NP plots.
        new_fp.write(f'coastlines={coastlines} #1:pscoast, 2:pscoast+InSAR.\n')

    new_fp.write("#coast_color is gray82 for off-white, or gray10 for dark coastlines.\n")

    # Adjust max/min latitudes for mapping points, otherwise points on edge aren't visible.
    if results['options']['output_choice'] in [1, 4]:
        buffer = 5
        if results['maxlat'] <= 90 - buffer:
            results['maxlat'] += buffer
        if results['minlat'] >= -90 + buffer:
            results['minlat'] -= buffer

    # Record text formats, which are the same for all projections and data types.
    new_fp.write('title_format="0 0 30 0 0 MC"\n')
    new_fp.write('blurb_format="0 0 15 0 1 ML"\n')
    new_fp.write('units_format="0 0 13 0 0 MC"\n')
    # Record units for the scale, which are the same for all projections and data types.
    new_fp.write(f"scale_units=\"{results['units'][i]}\"\n")

    new_fp.write('misc_range=" -R0/1/0/1 -JX1c "\n')
    new_fp.write("#grdcut requires actual limits, but if grdimage uses them: GMT Fatal Error: grdimage could not allocate memory [21.69 Gb, n_items = 5823567396]\n")
    new_fp.write('minlon=%.3f\n' % 0.0)  # results['minlon']
    new_fp.write('maxlon=%.3f\n' % 360.0)  # results['maxlon']
    new_fp.write('minlat=%.3f\n' % -90.0)  # results['minlat']
    new_fp.write('maxlat=%.3f\n' % 90.0)  # results['maxlat']

    if title:
        # All projections are always available.
        if i == 0:  # Only print this guide for the first map.
            new_fp.write("#Global projections:\n")
            new_fp.write("#    1 - Robinson\n")
            new_fp.write("#    2 - Winkel Tripel\n")
            new_fp.write("#    3 - Mollweide\n")
            new_fp.write("#    4 - Miller\n")
            new_fp.write("#Polar projections:\n")
            new_fp.write("#  101 - N. Azimuthal Equidistant\n")
            new_fp.write("#  102 - S. Azimuthal Equidistant\n")
            new_fp.write("#Specific regions:\n")
            new_fp.write("# 1001 - North America\n")
            new_fp.write("# 1002 - South America\n")
            new_fp.write("# 1003 - Africa\n")
            new_fp.write("# 1004 - Greenland\n")
            new_fp.write("# 1005 - South Asia\n")
            new_fp.write("# 1006 - Australia\n")
            new_fp.write("# 1007 - Europe\n")
            new_fp.write("# 1101 - Contiguous United States\n")
            new_fp.write("# 1102 - California\n")

        if polar == 0:
            new_fp.write(f"{'#' if i > 0 else ''}projection_choice={plot_options.projection}\n")
        else:
            if polar == 1:  # North pole.
                new_fp.write(f"{'#' if i > 0 else ''}projection_choice=101\n")
            else:  # South pole.
                new_fp.write(f"{'#' if i > 0 else ''}projection_choice=102\n")

        new_fp.write("standard_circle=0 #1=all specific regions use standard circular projection.\n")
        new_fp.write("standard_rect=0 #1=all specific regions use standard rectangular projection.\n")
        # . ./projections.sh -> here you would need to execute your shell script
        new_fp.write("if [ $projection_choice == 101 ]\n")
        new_fp.write("then\n")

        minlat = results['minlat'] if polar == 1 else 0.0
        new_fp.write(f"  minlat={minlat:.3f}\n")
        new_fp.write("  actual_range=\" -R0.0/360.0/$minlat/90.0 \"\n")
        polar_radius = 90 - minlat
        new_fp.write(f"  polar_radius={polar_radius}\n")
        new_fp.write(f"  projection=\" -JE0/90.0/{polar_radius}/${{map_width}}c \" #N. Azimuthal Equidistant\n")
        new_fp.write("elif [ $projection_choice == 102 ]\n")
        new_fp.write("then\n")

        maxlat = results['maxlat'] if polar == 2 else 0.0
        new_fp.write(f"  maxlat={maxlat:.3f}\n")
        new_fp.write("  actual_range=\" -R0.0/360.0/-90.0/$maxlat \"\n")
        polar_radius = 90 + maxlat
        new_fp.write(f"  polar_radius={polar_radius}\n")
        new_fp.write(f"  projection=\" -JE0/-90.0/{polar_radius}/${{map_width}}c \" #S. Azimuthal Equidistant\n")
        new_fp.write("fi\n")

        new_fp.write(f"range=\" -R${{minlon}}/${{maxlon}}/${{minlat}}/${{maxlat}} \"\n")
        new_fp.write(f"map_pos=\" -Xa${{map_x}}c -Ya${{map_y}}c \"\n")

        if len(results['rgb']) == 3 and results['rgb_choice'] >= 2:
            new_fp.write("scale_width=1.2 #Override for RGB maps.\n")

        new_fp.write("scale_pos=\" -D${{scale_x}}c/${{scale_y}}c/${{scale_length}}c/${{scale_width}}c \"\n")
        new_fp.write("units_x=$(bc <<< \"scale=5; $scale_x+$scale_width/2\")\n");
        new_fp.write("units_y=$(bc <<< \"scale=5; $scale_y+$scale_length/2\")\n");
        new_fp.write("units_pos=\" -Xa${{units_x}}c -Ya${{units_y}}c \"\n")
        new_fp.write("blurb_pos=\" -Xa${{blurb_x}}c -Ya${{blurbs_y}}c \"\n")
        new_fp.write("blurb2_pos=\" -Xa${{blurb2_x}}c -Ya${{blurbs_y}}c \"\n")
    else: print("!!!WARNING!!! NO TITLE!")

    # Plot data, with title on top.
    new_fp.write(f"title=\"{title}\"\n")
    if len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
        new_fp.write(f"$gmt_prefix grdimage red.nc green.nc blue.nc $boundary $resolution $range $projection $map_pos $start > $plot_base.ps\n")
    else:
        new_fp.write(f"$gmt_prefix grdimage $data_name $boundary $resolution $range $projection $map_pos -Cmap.cpt $start > $plot_base.ps\n")

    write_gmt_coastlines(new_fp)

    # If marker_lats multivectors have the right size, plot markers in the color "sandy brown".
    # Plot markers before mascons because mascons are smaller than markers.
    if len(results['marker_lats']) == len(results['latlon']['outputs']) and len(results['marker_lons']) == len(results['latlon']['outputs']) and results['latlon']['outputs']:
        if results['marker_lats'][i] and results['marker_lons'][i]:
            new_fp.write(f"$gmt_prefix psxy -N $data_name -bcmarker_lons/marker_lats -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n")

    new_fp.write("#Uncomment to put a marker at echoed coords, given as lon lat:\n")
    new_fp.write("#echo -85.19 -77.36 | $gmt_prefix psxy -N -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n")

    # If requested and mascon_lats/lons vectors aren't empty, plot mascon centers in the color "saddle brown".
    if plot_options['plot_mascons'] != 0 and results['latlon']['mascon_lats']:
        new_fp.write(f"$gmt_prefix psxy $data_name -bcmascon_lons/mascon_lats -Sc0.01c -G139/69/19 $range $projection $map_pos $middle >> $plot_base.ps\n")

    new_fp.write("if [ $montage != 0 ]\n")
    new_fp.write("then\n")
    new_fp.write("  title=${prefixes[$index]}\" \"$title\n")
    new_fp.write("  title_format=\"0 0 30 0 0 ML\" #Left-justify so montage titles are uniform.\n")
    new_fp.write("  title_x=$(bc <<< \"scale=5; $blurb_x-0.1\")\n")
    new_fp.write("else\n")
    new_fp.write("  title_x=$(bc <<< \"scale=5; $map_x+$map_width/2\")\n")
    new_fp.write("fi\n")
    new_fp.write("title_pos=\" -Xa${title_x}c -Ya${title_y}c \"\n")
    new_fp.write("echo $title_format $title | $gmt_prefix pstext -N $title_pos $misc_range $middle >> $plot_base.ps\n")

    # Draw color scale with units printed above.
    write_gmt_colorscale(new_fp, results, 0)  # 0 = next to map for GMT PDF output.

    # Print blurb about data range, or masked amplitudes for phase plots.
    # Don't print blurb if disabled or if this is an RGB map.
    if len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
        plot_options['blurb_disabled'] = 1
    if plot_options['blurb_disabled']:
        new_fp.write("blurb_contents=\"\"\n")

    new_fp.write("echo $blurb_format $blurb_contents | $gmt_prefix pstext -N $blurb_pos $misc_range $middle >> $plot_base.ps\n")

    # Only print error bars if they're the right length, and this one is > 0.0.
    # Need to have separate copies for unstructured and latlon maps because
    # in each case the "right length" is defined by a different multivector.
    blurb2_written = 0
    if len(results['latlon']['outputs']) == len(results['error_bars']) and results['latlon']['outputs']:
        if results['error_bars'][i] > 0.0:
            blurb2_written = 1
            new_fp.write(f"blurb2_contents=\"Error bar: {results['error_bars'][i]:.1f} $scale_units\"\n")

    if not blurb2_written:
        new_fp.write("blurb2_contents=\"\" #Error bar: N/A $scale_units\n")

    new_fp.write("echo $blurb_format $blurb2_contents | $gmt_prefix pstext -N $blurb2_pos $misc_range $end >> $plot_base.ps\n")





def write_gmt_defs(new_fp):
    """
    Purpose: This function writes some clarifying definitions that help
            me to consistently write working GMT data plotting commands.
            2017-04-04 Update: Also added GMT4/5 compatibility which uses
                            the "GMT5" define from definitions.hpp.
                            THIS ISN'T FINISHED YET!
    Input:  new_fp - file handle to current gmt script file
    """
    new_fp.write("#######################################################\n")
    new_fp.write("#Clarifying definitions. Do not change!################\n")
    new_fp.write("start=\" -K \" #Should always redirect using > to write new PS.\n")
    new_fp.write("middle=\" -O -K \" #Should always redirect using >> to append to PS.\n")
    new_fp.write("end=\" -O \" #Should always redirect using >> to append to PS.\n")
    new_fp.write("#######################################################\n")

    # Assuming GMT5 as a global variable here.
    if plot_options['GMT5']:
        new_fp.write("gmt5=1 #1/0 = GMT v5/4. GMTv5 support not finished.\n")
    else:
        new_fp.write("gmt5=0 #1/0 = GMT v5/4. GMTv5 support not finished.\n")

    new_fp.write("if [ $gmt5 != 0 ]\n")
    new_fp.write("then\n")
    new_fp.write("  gmt_prefix=\"gmt \"\n")
    new_fp.write("else\n")
    new_fp.write("  gmt_prefix=\"\"\n")
    new_fp.write("fi\n")









def finish_flip_backgrounds(flip_fp):
    """
    Purpose: This function finishes writing the flip_backgrounds.sh script.

    Input:  
        flip_fp - file pointer to flip_backgrounds.sh script file
    """
    flip_fp.write("for current_base in $all_bases\n")
    flip_fp.write("do\n")
    flip_fp.write("  #Change black to cyan temporarily.\n")
    flip_fp.write("  convert $current_base.png -fill cyan -opaque black $current_base.png\n")
    flip_fp.write("  #Change white to black.\n")
    flip_fp.write("  convert $current_base.png -fill black -opaque white $current_base.png\n")
    flip_fp.write("  #Change temporary cyan to white.\n")
    flip_fp.write("  convert $current_base.png -fill white -opaque cyan $current_base.png\n")
    flip_fp.write("done\n")

def finish_trim(trim_fp):
    """
    Purpose: This function finishes writing the trim.sh script.

    Input:  
        trim_fp - file pointer to trim.sh script file
    """
    trim_fp.write("for current_base in $all_bases\n")
    trim_fp.write("do\n")
    trim_fp.write("  convert $current_base.png -trim $current_base.png\n")
    trim_fp.write("done\n")

if __name__ == "__main__":
    write_gmt_scripts(plot_options, grid, results)