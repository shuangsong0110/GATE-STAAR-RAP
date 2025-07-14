# GATE-STAAR-RAP

This is the source code for the GATE-STAAR app that runs on the DNAnexus Platform. For more information about how to run or modify it, see https://documentation.dnanexus.com.

See also [GATE-STAAR](https://github.com/shuangsong0110/GATE-STAAR)

## Tutorial 1: Run GATE-STAAR in coding regions

### Step 1: Cloning an Applet
```
wget ...
```
### Step 2: Modify the R script file
Open the file `/gate-staar-coding/resources/home/dnanexus/STAAR_coding_backup.R`, then add your own DNAnexus token, and paths to the null model, varRatio, annotation catalog, and agds files.

### Step 3: Build the Applet
```
cd ./gate-staar-coding
dx build --overwrite
```
### Step 4: Run the Applet
```
for SLURM_ARRAY_TASK_ID in {1..22}
do
dx run gate_staar_coding -iarrayid=${SLURM_ARRAY_TASK_ID} --priority=low --instance-type="mem3_ssd1_v2_x4"
done
```


## Tutorial 2: Run GATE-STAAR in noncoding regions

### Step 1: Cloning an Applet
```
wget ...
```
### Step 2: Modify the R script file
Open the file `/gate-staar-noncoding/resources/home/dnanexus/STAAR_coding_backup.R`, then add your own DNAnexus token, and paths to the null model, varRatio, annotation catalog, and agds files.

### Step 3: Build the Applet
```
cd ./gate-staar-noncoding
dx build --overwrite
```
### Step 4: Run the Applet
```
for SLURM_ARRAY_TASK_ID in {1..22}
do
dx run gate_staar_noncoding -iarrayid=${SLURM_ARRAY_TASK_ID} --priority=low --instance-type="mem3_ssd1_v2_x4"
done
```
