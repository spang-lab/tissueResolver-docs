# In this script we run CIBERSORT in HiRes mode to determine cell type-specific proportions and gene expression
# for varying log foldchanges in our simulation

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
simulation="tissueResolver_${suffix}_logfc"
# vary logfc and percentage of modified genes fixed
ng=17
logfcs=( 0.5 0.8 1 1.5 2 )

# Step 0: Get the input files from the remote
mkdir $input_base_dir/$simulation
for logfc in "${logfcs[@]}"
do
    mkdir $input_base_dir/$simulation/$logfc
    for irun in $(seq 1 $number_runs)
    do
        mkdir $input_base_dir/$simulation/$logfc/$irun
        scp $server_name:$server_storage/bulks_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}.txt $input_base_dir/$simulation/$logfc/$irun
        scp $server_name:$server_storage/refsample_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}.txt $input_base_dir/$simulation/$logfc/$irun
    done
done

# # Step 1: Let CIBERSORT automatically build a signature from sc data for each logfc
# this we only need for the first run, because we use the same signature also for the other runs
for logfc in "${logfcs[@]}"
do
    irun=1
    input_dir=$input_base_dir/$simulation/$logfc/$irun
    echo $input_dir
    output_dir=$input_dir
    # run CIBERSORTin fractions mode using S-mode batch correction
    docker run \
        -v $input_dir:/src/data \
        -v $output_dir:/src/outdir \
        cibersortx/fractions \
        --username $mail \
        --token $token \
        --mixture bulks_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}.txt \
        --refsample refsample_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}.txt \
        --rmbatchSmode TRUE \
        --single_cell TRUE
done

for logfc in "${logfcs[@]}"
do
    irun=1
    output_dir=$input_base_dir/$simulation/$logfc/$irun
    scp $output_dir/CIBERSORTx_Adjusted.txt r4:$server_storage_cibersort/HiRes_fc_${logfc}_run_${irun}_${suffix}
done

# Copy the adjusted signature matrix to the other runs
for logfc in "${logfcs[@]}"
do
    for irun in $(seq 2 $number_runs)
    do
        input_dir=$input_base_dir/$simulation/$logfc/$irun
        output_dir=$input_base_dir/$simulation/$logfc/1
        cp $output_dir/CIBERSORTx_sigmatrix_Adjusted.txt $input_dir
    done
done

# Step 2: Determine cell type specific expression using the sigmatrix determined above
for logfc in "${logfcs[@]}"
do
    for irun in $(seq 1 $number_runs)
    do
        input_dir=$input_base_dir/$simulation/$logfc/$irun
        output_dir=$input_dir
        # run CIBERSORTin HiRes mode using B-mode batch correction
        docker run \
            -v $input_dir:/src/data \
            -v $output_dir:/src/outdir \
            cibersortx/hires \
            --username $mail \
            --token $token \
            --mixture bulks_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}.txt \
            --sigmatrix CIBERSORTx_sigmatrix_Adjusted.txt \
            --rmbatchBmode TRUE
    done
done

# Step 3: copy results to the server for downstream analysis
for logfc in "${logfcs[@]}"
do
    for irun in $(seq 1 $number_runs)
    do
        output_dir=$input_base_dir/$simulation/$logfc/$irun
        ssh $server_name "mkdir $server_storage_cibersort/HiRes_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}"
        scp $output_dir/CIBERSORTxGEP_NA_Fractions-Adjusted.txt $server_name:$server_storage_cibersort/HiRes_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}
        scp $output_dir/CIBERSORTxHiRes_NA_*.txt $server_name:$server_storage_cibersort/HiRes_fc_${logfc}_pct_${ng}_run_${irun}_${suffix}
    done
done



