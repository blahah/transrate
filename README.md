This is a forked branch of Transrate 1.0.4 that should work with the newest versions of Salmon (v1.7.0 , 01/03/2022) and snap (v2.0.1, 01/03/2022). 
***(Important note: I'm running the program on linux and was unsure regarding the macOS version of snap currently in use and as such there is no download link for the macOS version of snap on the deps.yaml file)***

This forked brach was created after issues in read mapping between RNA samples and a de novo built transcriptome built with those same samples. Mapping values on available recent versions of transrate would be much lower than expected and when compared to older versions of transrate. Also, potential brigdes would always be 0 when running the current version transrate.

With the changes made in this branch, transrate ran as expected and mapping values were comparable to that of older versions as well and potential bridges are calculate once again.

This is fork is mostly based on the v1.0.4 fork by [dfmoralesb](https://github.com/dfmoralesb/transrate) and builds upon changes already present there. Instalation instructions are also similar

To install this unofficial version of transrate v1.0.4.1:
```
$ git clone https://github.com/pmomadeira/transrate.git

$ cd transrate

$ gem build transrate.gemspec

$ sudo gem install transrate-1.0.3.gem
```

This installation of transrate should run with conda installations of salmon and snap, as long as you run it in the proper environment. If issues arise try the following [adapted from vhfsantos response](https://github.com/blahah/transrate/issues/202#issuecomment-457245624):
 ```
 #Create a new conda environment
 $conda create --name myTransrateEnv
 
 #Activate the environment
 $conda activate myTransrateEnv
 
 #Install salmon and snap in the selected environment
 $conda install -c bioconda salmon
 $conda install -c bioconda snap-aligner
 
 #List the path to both salmon and snap
 $whereis salmon
/home/user/miniconda3/envs/myTransrateEnv/bin/salmon
$whereis snap-aligner
/home/user/miniconda3/envs/myTransrateEnv/bin/snap-aligner

#Go into the bin directory of you transrate installation and create paths to conda's salmon and snap
$cd /path/to/transrate/bin/
$ln -s /home/user/miniconda3/envs/myTransrateEnv/bin/salmon
$ln -s /home/user/miniconda3/envs/myTransrateEnv/bin/snap-aligner

```
After this transrate should run without issues. Other than using different dependencies, this version should work like the original transrate, so follow general guidelines as instructed in [Transrate](http://hibberdlab.com/transrate/)
