if [[ $(bc <<< "$upper_limit != 0.0") == 1 ]]
then #Obtain floor of base-10 exponent.
  abs_upper_limit=$upper_limit
  if [[ $(bc <<< "$upper_limit < 0.0") == 1 ]]
  then
    abs_upper_limit=$(bc <<< "scale=5; $upper_limit*(-1)")
  fi
  upper_exp=$(bc -l <<< "scale=5; l($abs_upper_limit)/l(10)")
  printf -v upper_exp_int "%.0f" "$upper_exp"
  if [[ $(bc <<< "$upper_exp_int > $upper_exp") == 1 ]]
  then
    upper_exp_int=$(bc <<< "scale=5; $upper_exp_int-1")
  fi
  upper_exp=$upper_exp_int
else
  upper_exp=0
fi
if [[ $(bc <<< "$lower_limit != 0.0") == 1 ]]
then #Obtain floor of base-10 exponent.
  abs_lower_limit=$lower_limit
  if [[ $(bc <<< "$lower_limit < 0.0") == 1 ]]
  then
    abs_lower_limit=$(bc <<< "scale=5; $lower_limit*(-1)")
  fi
  lower_exp=$(bc -l <<< "scale=5; l($abs_lower_limit)/l(10)")
  printf -v lower_exp_int "%.0f" "$lower_exp"
  if [[ $(bc <<< "$lower_exp_int > $lower_exp") == 1 ]]
  then
    lower_exp_int=$(bc <<< "scale=5; $lower_exp_int-1")
  fi
  lower_exp=$lower_exp_int
else
  lower_exp=0
fi
#First, see if either exponent is > $digits+1.
if [[ $(bc <<< "$upper_exp > $digits+1") == 1 || $(bc <<< "$lower_exp > $digits+1") == 1 ]]
then
  notation="e" #If limits are too big, numbers appear in scientific notation.
  data_min_print=$data_min_e
  data_max_print=$data_max_e
fi
#Next, see if either exponent is < -$digits.
if [[ $(bc <<< "$upper_exp < -$digits") == 1 || $(bc <<< "$lower_exp < -$digits") == 1 ]]
then
  notation="e" #If limits are too small, numbers appear in scientific notation.
  data_min_print=$data_min_e
  data_max_print=$data_max_e
fi
scale_format=" --D_FORMAT=%.${digits}${notation} "
