top_overflow=$(bc <<< "$data_max_f > $upper_limit")
bottom_overflow=$(bc <<< "$data_min_f < $lower_limit")
if [[ $top_overflow == 1 && $bottom_overflow == 1 ]]
then
  overflow=" -E " #-E=both,-Ef=top,-Eb=bottom
elif [[ $top_overflow == 1 ]]
then
  overflow=" -Ef "
elif [[ $bottom_overflow == 1 ]]
then
  overflow=" -Eb "
else
  overflow=""
fi
