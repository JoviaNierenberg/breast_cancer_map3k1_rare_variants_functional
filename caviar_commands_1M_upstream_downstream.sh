#!/bin/bash
#
#PBS -N CAVIAR_UD1K
#PBS -S /bin/bash
#PBS -o /zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/Upstream_Downstream_1M/caviar_job.out
#PBS -e /zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/Upstream_Downstream_1M/caviar_job.err
#PBS -l nodes=1:ppn=1
#PBS -l vmem=20gb,mem=20gb
#PBS -m abe
#PBS -M Jovia.Nierenberg@ucsf.edu

# navigate to caviar directory
cd /zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/Upstream_Downstream_1M/

# function to run caviar commands
run_all_caviar_for_parameter () {
	for BLOCK in $(seq 1 $2); do
		echo "block $BLOCK"
		R_MATRIX="r_matrix_rsq_$1_block_$BLOCK.txt"
		Z_FILE="z_file_rsq_$1_block_$BLOCK.txt"
		OUT_NAME="caviar_out_rsq_$1_block_$BLOCK"
		OUT_NAME_MULT="caviar_out_mult_causal_rsq_$1_block_$BLOCK"
		if [[ -f $R_MATRIX && -f $Z_FILE ]]; then
			echo "files exist for block, running caviar:"
			CAVIAR -o $OUT_NAME -l $R_MATRIX -z $Z_FILE -c 1
			CAVIAR -o $OUT_NAME_MULT -l $R_MATRIX -z $Z_FILE
		else
			echo "no files in block"
		fi	done
}

# run caviar commands for rsq thresholds 0.2 and 0.4
run_all_caviar_for_parameter 04 226
run_all_caviar_for_parameter 02 191