# Workshop1
This workshop is catered to all levels of familiarity/expertise with HiC and computer programming, especially beginners. The workshop is in Stewart S3/3 June 20th from 1-5pm.

## Before the workshop:
You will need the following programs: Google chrome. Python, MATLAB.\
To use python I recommend anaconda https://www.anaconda.com/download. You can use spyder or vscode IDEs within Anaconda Navigator to start coding in python. 

### Install the following python packages by typing the following commands in the console:
```
pip install cooler
pip install cooltools
pip install numpy
pip install pandas
pip install matplotlib
```
### Request access to dataset hosted by Zenodo.
Create a zenodo account: https://zenodo.org/ you can login using your github account.\
Confirm your email attached to your zenodo account by clicking the link in the email sent to you by zenodo.\
Go to: https://zenodo.org/record/8040563. \
Request access to the data.\
Once approved to access the dataset, choose your favourite organism and download the appropriate .zip file.\
Open the folder. You should see a .hic file in the folder.\
If your folder doesn't have a .hic file don't worry, download another organism from zenodo.

## Open https://aidenlab.org/juicebox/ in chrome.
Open your .hic file using juicebox. You should see a plot of your HiC matrix. You can use the dropdown menus to display data for different chromosomes.

### Install MATLAB

### Go to: https://csynth.github.io/csynth/csynth.html
Command + option + I (on mac/chrome) to open developer tools\
Control + shift + C (on windows/chrome) to open developer tools\
In the console type: gl.getParameter(gl.MAX_TEXTURE_SIZE)\
Remember this number.

!!!\
There is a bug in other instances of CSynth limiting the number of monomers in the 3D model < 999 if the .txt upload format is used.\
https://csynth.molbiol.ox.ac.uk/csynth/login \
!!!\

Use https://csynth.github.io/csynth/csynth.html instead.\

### Email me the number, whether you use windows or mac, and the organism you chose. My email is my first name dot my last name @mail.mcgill.ca

# Part 1 (30mins)
Introduction to HiC

### Clone this github repo.
# Part 2 (2hrs)
In this part of the workshop we will compute contact probability or P(s) curves for many organisms across the tree of life.

Use the data here:
https://zenodo.org/record/8040563 and the following python script: plot_P(s)_workshop.py

Make sure you save the P(s) curve data. Name the .csv file according to your organism. If you finish repeat for a different organism. Look around the room to see the shapes of P(s) curves other people are calculating. How are they similar/different to your organism? Does your organisms exhibit chromosome to chromosome variation?

### MAC ONLY:
If you want to visualize the 3D conformation of chromosomes from your organisms open a new terminal window and type:

In terminal:
```
pip install cooler

cd "your organism folder"

cooler dump -t chroms "name_of_mcool_file"::/resolutions/5000
```
You should see an output similar to:
```
HiC_scaffold_1
HiC_scaffold_2
HiC_scaffold_3
etc…
```
These are the names of the chromosomes in the .mcool file. Change the bash code below: change the length of the for loop depending on how many chromosomes you have, change the file names, and change the .mcool file name. Then execute the command in terminal.

```
for i in $(seq 1 84);
do cooler dump -t pixels -r HiC_scaffold_${i} -r2 HiC_scaffold_${i} -o temp.txt --join GSM5182720_LetJap1.0_HiC.mcool::/resolutions/5000;
cut -f 2,5,7 temp.txt | column -t -s $"\t" > Arctic_lamprey_muscle_HiC_scaffold_${i}.txt;
rm -f temp.txt;
done
```
Drag and drop your .txt file into CSynth browser window: https://csynth.github.io/csynth/csynth.html

# Part 3 (2hrs)
In this part of the workshop we will optomize CSynth parameters to give us the most accurate 3D conformations possible. To do this we will use a data for IMR90 human immune cells. There are two types of data here: HiC data, and FISH (fluorescent in-situ hybridization) data. For these cells, the xyz position of many pairs of fluorescence probes that hybridize to a specific sequence of the genome, were recorded. We can compare the euclidian distance between fluorescence probes measured by microscopy to the euclidian distance between the same regions of the genome predicted by CSynth. We will change CSynth parameters to minimize these discrepancy as much as we can. We will repeat this analysis for three chromosomes: chr20, chr21, chr22.

### We will optomize two parameters that affect the 3D conformation significantly:
SPRINGPOW - this changes the relative importance of long and short distance effects. We will set SPRINGPOW to either (-2,-1,0) \
CONTACTFORCE - this changes magnitude of attraction between HiC contacts. We will set CONTACTFORCE to either (20,40,60,80,100)

Write your combination of parameters on the whiteboard in front of the room so it is not repeated by someone else.

## Here is a description of the files you will need:

The source publications i.e. where the data comes from:
DATA SOURCES.txt

This is the HiC data, unzip these files: \
chr20_5kb.RAWobserved.txt.zip \
chr21_5kb.RAWobserved.txt.zip \
chr22_5kb.RAWobserved.txt.zip

These are the primary sequence locations of the fluorescence probes: \
FISH_chr20bp.csv \
FISH_chr21bp.csv \
FISH_chr22bp.csv

These are the measured xyz positions of the fluorescence probes in the cell: \
FISH_chr20xyz.csv \
FISH_chr21xyz.csv \
FISH_chr22xyz.csv

These are the structures predicted by CSynth using the default parameters: \
IMR90_5kbp_chr20_MAPQG0_RAW_3D.xyz \
IMR90_5kbp_chr21_MAPQG0_RAW_3D.xyz \
IMR90_5kbp_chr22_MAPQG0_RAW_3D.xyz

Drag and drop your chr20_5kb.RAWobserved.txt into CSynth browser window: https://csynth.github.io/csynth/csynth.html /
Change SPRINGPOW and CONTACTFORCE parameters. Wait for the conformation to converge. For the CSynth structure to converge, the browser window needs to be displayed on your screen.

### Select:
Autoscale \
BED data source: rainbow \
Ribbon -> diameter: 30 \
Extras -> scripts -> eigen 

### Then Select:
save/load -> save_xyz
### Then select:
save/load -> save big image no gui

!!! save big image no gui removes the gui so impossible to use the gui after to click save_xyz !!!

Rename your structure/image according to the chromosome and parameters you used.

### Open and run: CSynth_FISH_Compare.m

Calculate the mean relative error (RE) between FISH probe separation and CSynth separation:
```
%RE = |H_ij - F_ij|/F_ij;
%F_ij FISH distance, H_ij HiC distance, ij genomic loci
```

Calculate the mean pearson correlation between FISH/HiC matricies.

Save this data in an excel spreadsheet.

# Part 4 (30min)
We will decide how to pool our data together.






