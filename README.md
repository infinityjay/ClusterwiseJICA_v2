# ClusterwiseJICA_v2

This project is the implementation for the Clusterwise J-ICA model including the simulation study scripts.

## New features

Based on the original Clusterwise J-ICA model, we mainly improve the fixed component number limitation.

* Varying the component number across clusters;

* Apply the VAF threshold rule and CHull to calculate the component number;

* Update the Tucker calculation method, for the situation when sources number larger than 8 apply a faster but not so precise calculation method.

## Structure of the project
As is shown in the following tree figure. There are 4 folders and each folder has specific functionality.

`analysis`: files to run the analysis based on the results from the simulation study. Both run on HPC or local.

`collect_result`: files to store the collected result from the result files of simulation study. Only run on HPC.

`functions`: functions for data simulation, Clusterwise J-ICA model and ILS based Clusterwise J-ICA model. The new Clusterwise J-ICA model with varying Q across clusters is implemented in the file *ClusterwiseJICA_Q.R* which is also the main model used in this thesis. The file *CJICA_asis_functions.R* includes almost all the functions used in the model and also the functions for ARI and Tucker calculation.

`scripts`: all the script files including the data generation, the job submission, the job analysis and the result check files. Can only be run on the HPC and the detailed information will be introduced in the following section.

```{bash}
.
├── LICENSE
├── README.md
├── analysis
│   ├── result_analysis.Rmd
├── collect_result
│   ├── comprehensive_results.RData
│   └── comprehensive_results.csv
├── functions
│   ├── CJICA_asist_functions.R
│   ├── ClusterwiseJICA.R
│   ├── ClusterwiseJICA_varyQ.R
│   ├── ILS_CJICA.R
│   ├── data_simulation_functions.R
│   └── icsfast.R
└── scripts
    ├── collect_result.R
    ├── dataGeneration.R
    ├── err_filter.sh
    ├── installPackages.R
    ├── run_analysis.R
    ├── submitJobs.sh
    └── testSimulationData.R

```


## Usage on the ALICE HPC
I use the GitHub to manage the whole project, and when I log into the ALICE platform, I also use `git clone` to pull the project. To avoid the large files, I will ignore all the generated RData file, except for the collected result file.

In the following sections, I will describe the steps to run the simulation step by step.

### Connecting to the ALICE
I mainly use the vscode to manipulate the files on the server, and I also use vscode to connect to the remote HPC servers. There is also an introduction about how to connect to ALICE through vscode on the webpage of ALICE wiki: https://pubappslu.atlassian.net/wiki/spaces/HPCWIKI/pages/37028145/Setting+up+VSCode+to+work+on+the+cluster.

Here is my settings to connect to the ALICE:

1. Setup hosts of ssh
In the ssh config file (usually `~/.ssh/config` on mac), add the following informations:
```bash
Host alice
    User <username> # your username, usually the student number
    Hostname login1.alice.universiteitleiden.nl # the login point, can also be login2.
    ProxyJump <username>@sshgw.leidenuniv.nl:22 # the jummp server, can also choose other url as jump server
    ServerAliveInterval 60
```

2. Config vscode
In the plugin page of vscode, you can search the plugin `Remote-SSH` and install it. Click on the remote ssh connection plugin, it will read the ssh config file by default, and you can just click the server `alice` which is configured in the step 1. Following the instruction, you will login the jump server and ALICE server, then you can connect to the ALICE successfully.

### Data generation
We separate the data generation and result analysis into two steps. And we will need to generate the data first.

Run the command `./scripts/dataGeneration.R` directly, and the generated data will be stored into the folder `simulation_data`. You can change the manipulate parameters by yourself in the script.

### Main job
The main job includes the job submission and the job analysis. Before submitting the jobs, you should install the R packages first. To install the packages, you should first load the R module use the command `module load R/4.4.0-gfbf-2023a`, then run the installing script `Rscript ./scripts/installPackages.R`. All the required packages are in the file `installPackages.R`, you can add the packages by yourself if needed.

The job submission script is the file `submitJobs.sh`, you can run it by command `./scripts/submitJobs.sh`. All the configures of the computing nodes are in the file, you can adjust it accordingly.

The main analysis script is the file `run_analysis.R`. For this file, you do not need to run it separately, and it already is executed by the job submission script. You can also change the analysis process accordingly.

Right now, the ALICE platform only support max 192 CPU nodes used by one user, you can rewrite the parallel job script accordingly.

### Results processing
You can see all the result files in the folder `./results`, and all the logs and error logs in the folder `./logs` which are all ignored in the `.gitignore` file. To check if there is any error in the all log files, you can execute the file `./scripts/err_filter.sh `.

After the error check, you can collect the results into one .csv file or .RData file. To execute the 
result collection process, you can run the command `Rscript ./scripts/collect_result.R` and all the results will be collected into .csv file or .RData in the folder `./collect_result`. You can also define which kind of factors you want to include into the result file by rewrite the `collect_all_results` function.

Finally, you can analyze the result from the data in the result file. I stored my analysis scripts in the folder `./analysis`.
### Some useful commands

```bash
find ./results/ -type f | wc -l # check completed result files 

squeue --user=<username> | grep ' R ' | wc -l # check all the running jobs for user xxx

scancel --user=<username> # cancel all the jobs submitted by the user

module load R/4.4.0-gfbf-2023a # load the R module for running R command

```
