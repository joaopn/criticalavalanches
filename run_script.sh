#!/bin/bash


if [ "$1" = "1" ]

then

	python3 ana/analyze_sim_ga.py --state subcritical -b 2 --de 4 --reps 50 --datafolder dat/NST/gamma_orlandi_01/dat/orlandi/ --ga 1,1.5,2

	python3 ana/analyze_sim_ga.py --state subcritical -b 2 --de 4 --reps 50 --datafolder dat/NST/gamma_orlandi_01/dat/gauss/ --ga 1,1.5,2


elif [ "$1" = "2" ]

then

	python3 ana/analyze_sim_ga.py --state reverberating -b 2 --de 4 --reps 50 --datafolder dat/NST/gamma_orlandi_01/dat/orlandi/ --ga 1,1.5,2

	python3 ana/analyze_sim_ga.py --state reverberating -b 2 --de 4 --reps 50 --datafolder dat/NST/gamma_orlandi_01/dat/gauss/ --ga 1,1.5,2

elif [ "$1" = "3" ]

then

	python3 ana/analyze_sim_ga.py --state critical -b 2 --de 4 --reps 50 --datafolder dat/NST/gamma_orlandi_02/dat/orlandi/ --ga 1,1.5,2

	python3 ana/analyze_sim_ga.py --state critical -b 2 --de 4 --reps 50 --datafolder dat/NST/gamma_orlandi_02/dat/gauss/ --ga 1,1.5,2


elif [ "$1" = "4" ]

then

	python3 ana/analyze_sim_ga.py --state subcritical -b 2 --de 4 --reps 50 --datafolder dat/NST/random_top_02/dat/random_local/ --ga 1,1.5,2

	python3 ana/analyze_sim_ga.py --state subcritical -b 2 --de 4 --reps 50 --datafolder dat/NST/random_top_02/dat/random_nonlocal/ --ga 1,1.5,2


elif [ "$1" = "5" ]

then

	python3 ana/analyze_sim_ga.py --state reverberating -b 2 --de 4 --reps 50 --datafolder dat/NST/random_top_02/dat/random_local/ --ga 1,1.5,2

	python3 ana/analyze_sim_ga.py --state reverberating -b 2 --de 4 --reps 50 --datafolder dat/NST/random_top_02/dat/random_nonlocal/ --ga 1,1.5,2


elif [ "$1" = "6" ]

then

	python3 ana/analyze_sim_ga.py --state critical -b 2 --de 4 --reps 50 --datafolder dat/NST/random_top_02/dat/random_local/ --ga 1,1.5,2

	python3 ana/analyze_sim_ga.py --state critical -b 2 --de 4 --reps 50 --datafolder dat/NST/random_top_02/dat/random_nonlocal/ --ga 1,1.5,2

fi
