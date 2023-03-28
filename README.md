# Building Moore

To build Moore, follow the steps [here](https://gitlab.cern.ch/rmatev/lb-stack-setup) up through compile. Use the latest version of the master branch. If you would like to use a nonstandard branch, follow instructions [here](https://gitlab.cern.ch/lhcb/Moore/-/blob/bnoc_run3/Hlt/Hlt2Conf/python/Hlt2Conf/lines/bnoc/README.md) but be aware that it may not run--for example, the `bnoc_run3` branch is often out of date. To update to master, I used the command 

```
git merge origin master
```
Best practice would be to develop with a local copy of the master branch. New code should now be pushed to `master`, whereas before it was pushed to `bnoc_run3`.

# Writing the HLT2 Line

The run 2 HLT lines live [here](https://gitlab.cern.ch/lhcb/Hlt/-/tree/2018-patches/Hlt/Hlt2Lines/python/Hlt2Lines), and the B-->K pi0 line lives [here](https://gitlab.cern.ch/lhcb/Hlt/-/tree/2018-patches/Hlt/Hlt2Lines/python/Hlt2Lines/B2Kpi0). New lines go in the Moore project with path
```
Moore/Hlt/Hlt2Conf/python/Hlt2Conf/lines
```
and should be booked in 

```
Hlt/Hlt2Conf/python/Hlt2Conf/lines/bnoc/hlt2_bnoc.py
```
Note this is not necessary to run the code, it only books it with `@register_line`. You will notice when writing the line that the dependencies are very different. What I did was I followed the [tutorial](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/hlt2_line.html#) to write an Hlt2 line and modified it based on the run 2 line. There is a B2OC which is fully implemented with ThOr functors [here](https://gitlab.cern.ch/lhcb/Moore/-/tree/master/Hlt/Hlt2Conf/python/Hlt2Conf/lines/b_to_open_charm_thor) which is helpful as an example. 

One obvious difference is that we need ThOr functors instead of LoKi, this conversion is described [here](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/thor_transition.html) and the comprehensive list of ThOr functors is [here](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/selection/thor_functors_reference.html). If any ThOr functors are missing, you can open an issue [here](https://gitlab.cern.ch/lhcb-dpa/project/-/issues/61). 

Make sure to use ```make_pvs_v2()``` rather than ```make_pvs()```, which will throw an error. If one or more of the decay products is neutral, when returning ```ParticleCombiner```, specify that ```ParticleCombiner="ParticleAdder"```. For example, mine reads
```
return ParticleCombiner(
        [kaons, pions],
        ParticleCombiner="ParticleAdder",
        DecayDescriptor="[B+ -> K+ pi0]cc",
        CombinationCut=combination_code,
        CompositeCut=composite_code,
    )
```

Once the line is complete, a branch and merge request should be opened in the relevant branch, following the guidelines [here](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/selection/hlt2_guidelines.html#summary). As an example, here is my 

# Running 
To run the code, you must write an options file. I put mine in the directory 'stack' which contains all the projects, but eventually it will go into a directory containing all the options files 
```
Hlt/Hlt2Conf/options
```
Then, follow this [tutorial](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/hlt2_line.html#running) which will tell you how to write an options file that runs over the minimum bias input. 

To run over MC files, find the file you want in bookkeeping by following [this tutorial](https://lhcb.github.io/starterkit-lessons/first-analysis-steps/bookkeeping.html), and downloading it as a python file. Using the path to that file, create an array in your options file which contains the paths to the files you want. The path should be 
```
root://eoslhcb.cern.ch//eos/lhcb/grid/prod/INSERT PATH YOU COPIED FROM BOOKKEEPING HERE
```
You may encounter MC files of different formats, [this tutorial](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/different_samples.html) will show you how to handle them--but it won't tell you everything! Also, some of the example options files it gives are inconsistent with what the tutorial says--for example options.input_raw_format is set in the tutorial but not in the example. I used the tutorial setting. It is also missing some things. You need to include two options not mentioned in the tutorial which are ```options.dddb_tag``` and ```options.conddb_tag```. The settings for these are specific to each decay, you can find them by following the instructions [here](https://lhcb.github.io/starterkit-lessons/first-analysis-steps/minimal-dv-job.html) under the heading "Database tags". Errors related to the decoding version can be handled following the instructions [here](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/different_samples.html#ft-decoding-version).

Most recently (at time of writing) it appears that LFNs are not accepted as input files, and they must be converted to PFNs. To convert to PFN, use

```
lb-dirac dirac-bookkeeping-genXMLCatalog --Options=LFN_from_bookkeeping.py --NewOptions=PFN.py
```
To run locally, use the command

```
./Moore/run gaudirun.py /path/to/options.py /path/to/PFN.py
```
 
# Options File

To run over more than one file, there are a few options. At first I was running over files in eos, but I realized not all the MC files were in there. Then, I looked into using a catalog instead, follwing [these](https://lhcb.github.io/starterkit-lessons/first-analysis-steps/files-from-grid.html) instructions under the "read files remotely..." heading. An example of using inputs from eos is [here](https://github.com/umd-lhcb/hlt-umd/blob/139e95f6c6d74769452e5a81d88a84fc5feaf269/test_b_to_kpi_xdigi.py).

# Moore and DaVinci

I only ran DaVinci on lxplus since the LDST file produced after applying my HLT only had a few hundred events. I created the DaVinci script following [this example](https://gitlab.cern.ch/lhcb/DaVinci/-/blob/master/DaVinciExamples/python/DaVinciExamples/tupling/option_davinci_configFuntuple.py). There are many examples floating around out there, the tutorials are very out of date on the subject of run 3, and there are lots of defunct naming conventions too. The example I used is working at the time of writing this, and my own DaVinci script is named 
option_davinci_tupling_from_hlt2.py. To run the code, I used the command
```
DaVinci/run davinci run-mc --inputfiledb Hlt2Output DV_scripts/data.yaml --joboptfile DV_scripts/option_davinci_hlt2.yaml --user_algorithms DV_scripts/option_davinci_tupling_from_hlt2:main
```
`data.yaml` specifies the inputs, `option_davinci_tupling_from_hlt2` is the python file I just described, and `option_davinci_hlt2.yaml` gives some output and other options. An example lives in this directory as well as the one mentioned earlier which I used to write my own.

Notably, at the time of writing this document you could not put any functor which uses the distance calculator as a variable in the ntuples, according to [this](https://mattermost.web.cern.ch/lhcb/pl/opuwpjo9qid3j8nx1qn8ts9ewe) discussion. This is being fixed as I write. 

# Moore and Ganga

I have a modified version of Moore and Rec that I need to use when submitting jobs to Ganga. For this, you must create the nightly build of Moore following instructions [here](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/ganga.html#build). Mine requires a modified version of Rec as well, so inside the build of Moore, I used the commands

```
git lb-use Rec
git lb-checkout Rec/Master Phys
```
For more information about the `lb-dev` environment, you can look [here](https://lhcb-core-doc.web.cern.ch/lhcb-core-doc/GitForLHCbUsers.html#using-git-for-lhcb-development). The higher level project should be built using the nightly, and any dependencies can be used and checked out inside that directory. From the nightly directory, you can change whatever you need and submit a job in Moore to ganga. [These instructions](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/ganga.html) describe how to write the script for Ganga, and [this](https://lhcb.github.io/starterkit-lessons/first-analysis-steps/davinci-grid.html) provides more information on writing the script and setting up the job. Overall, the starterkit has some documentation on Ganga. In this directory I have `MC_job.py` and `minbias_job.py` which provide examples of how to write Ganga scripts for Moore jobs. 

The file you need to specify inputs is the one you download by following the instructions from bookkeeping [here](https://lhcb.github.io/starterkit-lessons/first-analysis-steps/bookkeeping.html).

You need both the output LDST file as well as JSON file created along with the LDST in order to run DaVinci, instructions for how to produce this JSON file are [here](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/hlt2_analysis.html).

Do not try to use `lb-...` commands for pushing to a branch. They are defunct. For development use `git` commands from the project you are developing (for example, from `Rec` if I am doing functor development). The only purpose of these `lb-dev` environments as far as I can tell is to provide a copy of your code to Ganga so it can run a job on your version of Moore. 

As of today (3/28/23) the working Ganga script is titled `minbias_job.py`, and the compatible options file is titled `test_b_to_kpi.py`. To run the job all you need to do is 

```
lhcb-proxy-init
ganga minbias_job.py
```