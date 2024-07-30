# In this script we run CIBERSORT in fractions mode to determine cell type-specific proportions
# for simulated bulk data

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

# name of the file storing the proportions output by CIBERSORT
origin_file_name=CIBERSORTx_Adjusted.txt
# name of the simulation
simulation=steen

# Let CIBERSORT automatically build a signature matrix from sc data
number_runs=5
algo=sc
output_dir=$output_base_dir/$simulation/$algo
input_dir=$input_base_dir/$simulation
mkdir -p $output_dir
for irun in $(seq 1 $number_runs)
do
    # run CIBERSORTin fractions mode using S-mode batch correction
    docker run \
        -v $input_dir:/src/data \
        -v $output_dir:/src/outdir \
        cibersortx/fractions \
        --username $mail \
        --token $token \
        --single_cell TRUE \
        --refsample refsample_run_${irun}.txt \
        --mixture bulks_run_${irun}.txt \
        --fraction 0 \
        --rmbatchSmode TRUE
    destination_file_name="CIBERSORTx_Adjusted_${algo}_run_${irun}.txt"
    # copy the CIBERSORT reult to the server for further downstream
    scp $output_dir/$origin_file_name $server_name:$server_storage/$destination_file_name
done
