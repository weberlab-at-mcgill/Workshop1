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
Go to: https://zenodo.org/record/8040563.\
Request access to the data.\
Once approved to access the dataset, choose your favourite organism and download the appropriate .zip file.\
Open the folder. You should see a .hic file in the folder.\
If your folder doesn't have a .hic file don't worry, download another organism from zenodo.

## Open https://aidenlab.org/juicebox/ in chrome.
Open your .hic file using juicebox. You should see a plot of your HiC matrix. You can use the dropdown menus to navigate to different chromosomes/HiC_SCAFFOLDS.

### Install MATLAB

### Go to: https://csynth.github.io/csynth/csynth.html
Command + option + I (on mac/chrome) to open developer tools\
Control + shift + C (on windows/chrome) to open developer tools\
In the console type: gl.getParameter(gl.MAX_TEXTURE_SIZE)\
Remember this number. \

!!!\
There is a bug in other instances of CSynth limiting the number of monomers in the 3D model < 999 if the .txt upload format is used.\
This includes:\
https://csynth.molbiol.ox.ac.uk/csynthstatic/latest/csynth.html/ \
as well as\
https://csynth.molbiol.ox.ac.uk/csynth/login \
!!!\
Use https://csynth.github.io/csynth/csynth.html instead.

### Email me the number, whether you use windows or mac, and the organism you chose. My email is my first name dot my last name @mail.mcgill.ca

# Part 1
Introduction to HiC

# Part 2
3D Visualization of HiC Data\

For the CSynth structure to converge, the browser window needs to be displayed on your screen. Ask me about code to periodicaly ping server to fix this and have CSynth run in the background.

BED data source: rainbow

Explain files uploaded to GitHub:\
chr20_5kb.RAWobserved\
chr21_5kb.RAWobserved\
chr22_5kb.RAWobserved\

chr20_5kb.RAWobserved.txt\
chr21_5kb.RAWobserved.txt\
chr22_5kb.RAWobserved.txt\

CSynth_FISH_Compare.m\
DATA SOURCES.txt\

FISH_chr20bp.csv\
FISH_chr21bp.csv\
FISH_chr22bp.csv\

FISH_chr20xyz.csv\
FISH_chr21xyz.csv\
FISH_chr22xyz.csv\

GSE63525_IMR90_README.rtf\
IMR90_5kbp_chr20_MAPQG0_RAW_3D.xyz\
IMR90_5kbp_chr21_MAPQG0_RAW_3D.xyz\
IMR90_5kbp_chr22_MAPQG0_RAW_3D.xyz\

In terminal:
pip install cooler

cd into your organism folder

cooler dump -t chroms  name_of_mcool_file::/resolutions/5000
HiC_scaffold_1
HiC_scaffold_2
HiC_scaffold_3
etc…

change the length of the for loop depending on how many chromosomes you have, the file names, and the .mcool file name

for i in $(seq 1 84);  do cooler dump -t pixels -r HiC_scaffold_${i} -r2 HiC_scaffold_${i} -o temp.txt --join GSM5182720_LetJap1.0_HiC.mcool::/resolutions/5000; cut -f 2,5,7 temp.txt | column -t -s $"\t" > Arctic_lamprey_muscle_HiC_scaffold_${i}.txt; rm -f temp.txt; done

automatic_CSynth.js\

# Part 3
P(s) curves
https://zenodo.org/record/8040563.\
plot_P(s)_workshop.py\







