#!/bin/bash
#Circle is default if both=0 so polar works, and if both=1 for consistency.
if [ $standard_rect == 0 ]
then #Standard circle.
  map_x=9.0
  map_y=2.5
  map_width=15.0
  title_y=19.3
  scale_y=10.0
  scale_length=12.0
  scale_width=0.6
  scale_x=$(bc <<< "scale=5; $map_x-2.0")
  blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
  blurb2_x=$(bc <<< "scale=5; $blurb_x+12.8")
  blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  boundary=" -B60g60/20g20"
  actual_range=" -R0.0/360.0/-90.0/90.0 "
fi
if [[ $standard_rect != 0 && $standard_circle == 0 ]]
then #Standard rectangle.
  map_x=5.6
  map_y=1.5
  map_width=22.0
  title_y=20.3
  scale_y=9.5
  scale_length=15.0
  scale_width=0.6
  scale_x=$(bc <<< "scale=5; $map_x-2.5")
  blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
  blurb2_x=$(bc <<< "scale=5; $blurb_x+21.5")
  blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  boundary=" -B45g45/30g30/" #Don't remember why this boundary ends in a "/"...
  actual_range=" -R0.0/360.0/-90.0/90.0 "
fi
if [ $projection_choice == 1 ]
then
  map_x=5.5
  map_y=4.5
  map_width=22.0
  title_y=18.0
  scale_y=10.1
  scale_length=9.5
  scale_width=0.6
  projection=" -JN0/${map_width}c " #Robinson
  scale_x=$(bc <<< "scale=5; $map_x-2.6")
  blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
  blurb2_x=$(bc <<< "scale=5; $blurb_x+23.7")
  blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  boundary=" -B60g60/20g20"
  actual_range=" -R0.0/360.0/-90.0/90.0 "
elif [ $projection_choice == 2 ]
then
  map_x=5.5
  map_y=3.5
  map_width=22.0
  title_y=18.6
  scale_y=10.21
  scale_length=11.0
  scale_width=0.6
  projection=" -JR0/${map_width}c " #Winkel Tripel
  scale_x=$(bc <<< "scale=5; $map_x-2.5")
  blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
  blurb2_x=$(bc <<< "scale=5; $blurb_x+22.4")
  blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  boundary=" -B60g60/20g20"
  actual_range=" -R0.0/360.0/-90.0/90.0 "
elif [ $projection_choice == 3 ]
then
  map_x=5.5
  map_y=5.0
  map_width=22.0
  title_y=18.0
  scale_y=10.5
  scale_length=9.0
  scale_width=0.6
  projection=" -JW0/${map_width}c " #Mollweide
  scale_x=$(bc <<< "scale=5; $map_x-2.5")
  blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
  blurb2_x=$(bc <<< "scale=5; $blurb_x+23.5")
  blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  boundary=" -B60g60/20g20"
  actual_range=" -R0.0/360.0/-90.0/90.0 "
elif [ $projection_choice == 4 ]
then
  map_x=5.6
  map_y=1.5
  map_width=22.0
  title_y=20.3
  scale_y=9.5
  scale_length=15.0
  scale_width=0.6
  projection=" -JJ0/${map_width}c " #Miller
  scale_x=$(bc <<< "scale=5; $map_x-2.5")
  blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
  blurb2_x=$(bc <<< "scale=5; $blurb_x+21.5")
  blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  boundary=" -B45g45/30g30/" #Don't remember why this boundary ends in a "/"...
  actual_range=" -R0.0/360.0/-90.0/90.0 "
elif [ $projection_choice == 1001 ] #North America
then
  minlon=190.0
  maxlon=310.0
  minlat=0.0
  maxlat=85.0
  actual_range=" -R$minlon/$maxlon/$minlat/$maxlat "
  midlon=$(bc <<< "scale=5; ($maxlon+$minlon)/2")
  midlat=$(bc <<< "scale=5; ($maxlat+$minlat)/2")
  if [ $standard_circle != 0 ] #Standard circular map over-ride.
  then
    projection=" -JG$midlon/$midlat/25.0/${map_width}c " #Orthographic
  elif [ $standard_rect != 0 ]
  then
    projection=" -JJ0/${map_width}c " #Miller
  else
    map_x=6.0
    map_y=3.5
    map_width=20.0
    title_y=19.3
    scale_y=10.0
    scale_length=12.0
    scale_width=0.6
    projection=" -JL$midlon/$midlat/$minlat/$maxlat/${map_width}c " #Lambert Conformal Conic
    scale_x=$(bc <<< "scale=5; $map_x-3.0")
    blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
    blurb2_x=$(bc <<< "scale=5; $blurb_x+20.0")
    blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  fi
  boundary=" -BWESn10g10/10g10"
  pscoast_res=" -N2 "$pscoast_res_orig
elif [ $projection_choice == 1002 ] #South America
then
  minlon=275.0
  maxlon=330.0
  minlat=-57.0
  maxlat=13.0
  actual_range=" -R$minlon/$maxlon/$minlat/$maxlat "
  midlon=$(bc <<< "scale=5; ($maxlon+$minlon)/2")
  midlat=$(bc <<< "scale=5; ($maxlat+$minlat)/2")
  if [ $standard_circle != 0 ] #Standard circular map over-ride.
  then
    projection=" -JG$midlon/$midlat/25.0/${map_width}c " #Orthographic
  elif [ $standard_rect != 0 ]
  then
    projection=" -JJ0/${map_width}c " #Miller
  else
    map_x=6.0
    map_y=3.5
    map_width=11.0
    title_y=19.3
    scale_y=10.0
    scale_length=12.0
    scale_width=0.6
    projection=" -JL$midlon/$midlat/$minlat/$maxlat/${map_width}c " #Lambert Conformal Conic
    scale_x=$(bc <<< "scale=5; $map_x-3.0")
    blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
    blurb2_x=$(bc <<< "scale=5; $blurb_x+20.0")
    blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  fi
  boundary=" -BWEsN10g10/10g10"
elif [ $projection_choice == 1003 ] #Africa
then
  minlon=-20.0
  maxlon=55.0
  minlat=-37.0
  maxlat=38.0
  actual_range=" -R$minlon/$maxlon/$minlat/$maxlat "
  midlon=$(bc <<< "scale=5; ($maxlon+$minlon)/2")
  midlat=$(bc <<< "scale=5; ($maxlat+$minlat)/2")
  if [ $standard_circle != 0 ] #Standard circular map over-ride.
  then
    projection=" -JG$midlon/$midlat/25.0/${map_width}c " #Orthographic
  elif [ $standard_rect != 0 ]
  then
    projection=" -JJ0/${map_width}c " #Miller
  else
    map_x=6.0
    map_y=3.5
    map_width=12.0
    title_y=19.3
    scale_y=10.0
    scale_length=12.0
    scale_width=0.6
    projection=" -JL$midlon/$midlat/$minlat/$maxlat/${map_width}c " #Lambert Conformal Conic
    scale_x=$(bc <<< "scale=5; $map_x-3.0")
    blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
    blurb2_x=$(bc <<< "scale=5; $blurb_x+20.0")
    blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  fi
  boundary=" -BWEsN10g10/10g10"
elif [ $projection_choice == 1004 ] #Greenland
then
  minlon=283.0
  maxlon=347.0
  minlat=59.0
  maxlat=85.0
  actual_range=" -R$minlon/$maxlon/$minlat/$maxlat "
  midlon=$(bc <<< "scale=5; ($maxlon+$minlon)/2")
  midlat=$(bc <<< "scale=5; ($maxlat+$minlat)/2")
  if [ $standard_circle != 0 ] #Standard circular map over-ride.
  then
    projection=" -JG$midlon/$midlat/25.0/${map_width}c " #Orthographic
  elif [ $standard_rect != 0 ]
  then
    projection=" -JJ0/${map_width}c " #Miller
  else
    map_x=6.0
    map_y=3.5
    map_width=16.0
    title_y=19.3
    scale_y=10.0
    scale_length=12.0
    scale_width=0.6
    projection=" -JL$midlon/$midlat/$minlat/$maxlat/${map_width}c " #Lambert Conformal Conic
    scale_x=$(bc <<< "scale=5; $map_x-3.0")
    blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
    blurb2_x=$(bc <<< "scale=5; $blurb_x+20.0")
    blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  fi
  boundary=" -BWESn10g10/10g10"
elif [ $projection_choice == 1005 ] #South Asia
then
  minlon=40.0
  maxlon=140.0
  minlat=0.0
  maxlat=50.0
  actual_range=" -R$minlon/$maxlon/$minlat/$maxlat "
  midlon=$(bc <<< "scale=5; ($maxlon+$minlon)/2")
  midlat=$(bc <<< "scale=5; ($maxlat+$minlat)/2")
  if [ $standard_circle != 0 ] #Standard circular map over-ride.
  then
    projection=" -JG$midlon/$midlat/25.0/${map_width}c " #Orthographic
  elif [ $standard_rect != 0 ]
  then
    projection=" -JJ0/${map_width}c " #Miller
  else
    map_x=6.0
    map_y=3.5
    map_width=22.0
    title_y=19.3
    scale_y=10.0
    scale_length=12.0
    scale_width=0.6
    projection=" -JL$midlon/$midlat/$minlat/$maxlat/${map_width}c " #Lambert Conformal Conic
    scale_x=$(bc <<< "scale=5; $map_x-3.0")
    blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
    blurb2_x=$(bc <<< "scale=5; $blurb_x+20.0")
    blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  fi
  boundary=" -B10g10/10g10"
elif [ $projection_choice == 1006 ] #Australia
then
  minlon=110.0
  maxlon=156.0
  minlat=-45.0
  maxlat=-5.0
  actual_range=" -R$minlon/$maxlon/$minlat/$maxlat "
  midlon=$(bc <<< "scale=5; ($maxlon+$minlon)/2")
  midlat=$(bc <<< "scale=5; ($maxlat+$minlat)/2")
  if [ $standard_circle != 0 ] #Standard circular map over-ride.
  then
    projection=" -JG$midlon/$midlat/25.0/${map_width}c " #Orthographic
  elif [ $standard_rect != 0 ]
  then
    projection=" -JJ0/${map_width}c " #Miller
  else
    map_x=6.0
    map_y=3.5
    map_width=15.0
    title_y=19.3
    scale_y=10.0
    scale_length=12.0
    scale_width=0.6
    projection=" -JL$midlon/$midlat/$minlat/$maxlat/${map_width}c " #Lambert Conformal Conic
    scale_x=$(bc <<< "scale=5; $map_x-3.0")
    blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
    blurb2_x=$(bc <<< "scale=5; $blurb_x+20.0")
    blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  fi
  boundary=" -B10g10/10g10"
elif [ $projection_choice == 1007 ] #Europe
then
  minlon=-28.0
  maxlon=40.0
  minlat=35.0
  maxlat=72.0
  actual_range=" -R$minlon/$maxlon/$minlat/$maxlat "
  midlon=$(bc <<< "scale=5; ($maxlon+$minlon)/2")
  midlat=$(bc <<< "scale=5; ($maxlat+$minlat)/2")
  if [ $standard_circle != 0 ] #Standard circular map over-ride.
  then
    projection=" -JG$midlon/$midlat/25.0/${map_width}c " #Orthographic
  elif [ $standard_rect != 0 ]
  then
    projection=" -JJ0/${map_width}c " #Miller
  else
    map_x=6.0
    map_y=3.5
    map_width=18.0
    title_y=19.3
    scale_y=10.0
    scale_length=12.0
    scale_width=0.6
    projection=" -JL$midlon/$midlat/$minlat/$maxlat/${map_width}c " #Lambert Conformal Conic
    scale_x=$(bc <<< "scale=5; $map_x-3.0")
    blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
    blurb2_x=$(bc <<< "scale=5; $blurb_x+20.0")
    blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  fi
  boundary=" -B10g10/10g10"
elif [ $projection_choice == 1101 ] #Contiguous United States
then
  minlon=230.0
  maxlon=300.0
  minlat=23.0
  maxlat=55.0
  actual_range=" -R$minlon/$maxlon/$minlat/$maxlat "
  midlon=$(bc <<< "scale=5; ($maxlon+$minlon)/2")
  midlat=$(bc <<< "scale=5; ($maxlat+$minlat)/2")
  if [ $standard_circle != 0 ] #Standard circular map over-ride.
  then
    projection=" -JG$midlon/$midlat/25.0/${map_width}c " #Orthographic
  elif [ $standard_rect != 0 ]
  then
    projection=" -JJ0/${map_width}c " #Miller
  else
    map_x=6.0
    map_y=3.5
    map_width=22.0
    title_y=19.3
    scale_y=10.0
    scale_length=12.0
    scale_width=0.6
    projection=" -JL$midlon/$midlat/$minlat/$maxlat/${map_width}c " #Lambert Conformal Conic
    scale_x=$(bc <<< "scale=5; $map_x-3.0")
    blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
    blurb2_x=$(bc <<< "scale=5; $blurb_x+20.0")
    blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  fi
  boundary=" -B10g10/10g10"
  pscoast_res=" -N2 "$pscoast_res_orig
elif [ $projection_choice == 1102 ] #California
then
  minlon=233.0
  maxlon=250.0
  minlat=30.0
  maxlat=45.0
  actual_range=" -R$minlon/$maxlon/$minlat/$maxlat "
  midlon=$(bc <<< "scale=5; ($maxlon+$minlon)/2")
  midlat=$(bc <<< "scale=5; ($maxlat+$minlat)/2")
  if [ $standard_circle != 0 ] #Standard circular map over-ride.
  then
    projection=" -JG$midlon/$midlat/25.0/${map_width}c " #Orthographic
  elif [ $standard_rect != 0 ]
  then
    projection=" -JJ0/${map_width}c " #Miller
  else
    map_x=10.5
    map_y=3.5
    map_width=12.0
    title_y=19.3
    scale_y=10.0
    scale_length=12.0
    scale_width=0.6
    projection=" -JL$midlon/$midlat/$minlat/$maxlat/${map_width}c " #Lambert Conformal Conic
    scale_x=$(bc <<< "scale=5; $map_x-3.0")
    blurb_x=$(bc <<< "scale=5; $scale_x-2.5")
    blurb2_x=$(bc <<< "scale=5; $blurb_x+13")
    blurbs_y=$(bc <<< "scale=5; $title_y-1.2")
  fi
  boundary=" -B10g10/10g10"
  pscoast_res=" -N2 "$pscoast_res_orig
fi
