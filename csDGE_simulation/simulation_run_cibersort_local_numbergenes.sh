# In this script we run CIBERSORT in HiRes mode to determine cell type-specific proportions and gene expression
# for varying number of modified genes in our simulation

# specifiy the home directory on your local machine
mac_home=
# specifiy the home directory on your remote machine
server_home=
# specify here the name of your remote server for scp results between local and remote storage
server_name=
# specify your mail adress registered at CIBERSORT
mail=
# specify your mail adress provided CIBERSORT
token=

# local input/output directory for cibersort files
input_base_dir=$mac_home
# remote storage to scp the results to
server_storage=$server_home
server_storage_cibersort=$server_home


number_runs=1
suffix="B"
# suffix="cd8"
# name of the output folder
simulation="tissueResolver_${suffix}_numbergenes"
# vary precentage of modified genes and leave logfc fixed
number_genes=( 2 5 10 17 20 30 40 50 60 70 80 90 95 )
logfc=0.8

# Step 0: Get the input files from the remote
mkdir $input_base_dir/$simulation
for ng in "${number_genes[@]}"
do
    mkdir $input_base_dir/$simulation/$ng
    for irun in $(seq 1 $number_runs)
    do
        mkdir $input_base_dir/$simulation/$ng/$irun
        scp $server_name:$server_storage/bulks_fc_${fc}_pct_${ng}_run_${irun}_${suffix}.txt $input_base_dir/$simulation/$ng/$irun
        scp $server_name:$server_storage/refsample_fc_${fc}_pct_${ng}_run_${irun}_${suffix}.txt $input_base_dir/$simulation/$ng/$irun
    done
done


# # Step 1: Let CIBERSORT automatically build a signature from sc data for each logfc
# this we only need for the first run, because we use the same signature also for the other runs
for ng in "${number_genes[@]}"
do
    irun=1
    input_dir=$input_base_dir/$simulation/$ng/$irun
    echo $input_dir
    output_dir=$input_dir
    # run CIBERSORTin fractions mode using S-mode batch correction
    docker run \
        -v $input_dir:/src/data \
        -v $output_dir:/src/outdir \
        cibersortx/fractions \
        --username $mail \
        --token $token \
        --mixture bulks_fc_${fc}_pct_${ng}_run_${irun}_${suffix}.txt \
        --refsample refsample_fc_${fc}_pct_${ng}_run_${irun}_${suffix}.txt \
        --rmbatchSmode TRUE \
        --single_cell TRUE
done

# Copy the adjusted signature matrix to the other runs
for ng in "${number_genes[@]}"
do
    for irun in $(seq 2 $number_runs)
    do
        input_dir=$input_base_dir/$simulation/$ng/$irun
        output_dir=$input_base_dir/$simulation/$ng/1
        cp $output_dir/CIBERSORTx_sigmatrix_Adjusted.txt $input_dir
    done
done

# Step 2: Determine cell type specific expression using the sigmatrix determined above
for ng in "${number_genes[@]}"
do
    for irun in $(seq 1 $number_runs)
    do
    input_dir=$input_base_dir/$simulation/$ng/$irun
    output_dir=$input_dir
    # run CIBERSORTin HiRes mode using B-mode batch correction
    docker run \
        -v $input_dir:/src/data \
        -v $output_dir:/src/outdir \
        cibersortx/hires \
        --username $mail \
        --token $token \
        --mixture bulks_fc_${fc}_pct_${ng}_run_${irun}_${suffix}.txt \
        --sigmatrix CIBERSORTx_sigmatrix_Adjusted.txt \
        --rmbatchBmode TRUE
    done
done

# Step 3: copy results to the server for downstream analysis
for ng in "${number_genes[@]}"
do
    for irun in $(seq 1 $number_runs)
    do
        output_dir=$input_base_dir/$simulation/$ng/$irun
        ssh $server_name "mkdir $server_storage_cibersort/HiRes_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}"
        scp $output_dir/CIBERSORTxGEP_NA_Fractions-Adjusted.txt $server_name:$server_storage_cibersort/HiRes_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}
        scp $output_dir/CIBERSORTxHiRes_NA_*.txt $server_name:$server_storage_cibersort/HiRes_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}
    done
done