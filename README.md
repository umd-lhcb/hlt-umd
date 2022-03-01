# Building Moore

To build Moore, follow the steps [here](https://gitlab.cern.ch/rmatev/lb-stack-setup) up through compile. Use the latest version of the master branch. If you would like to use a nonstandard branch, follow instructions [here](https://gitlab.cern.ch/lhcb/Moore/-/blob/bnoc_run3/Hlt/Hlt2Conf/python/Hlt2Conf/lines/bnoc/README.md) but be aware that there may be bugs--for example, I was not able to run the example line on the bnoc_run3 branch without updating to master. To update to master, I used the command 

```
git merge origin/master
```

# Writing the HLT2 Line

The run 2 HLT lines live [here](https://gitlab.cern.ch/lhcb/Hlt/-/tree/2018-patches/Hlt/Hlt2Lines/python/Hlt2Lines), and the B-->K pi0 line lives [here](https://gitlab.cern.ch/lhcb/Hlt/-/tree/2018-patches/Hlt/Hlt2Lines/python/Hlt2Lines/B2Kpi0). New lines go in the Moore project with path
```
Moore/Hlt/Hlt2Conf/python/Hlt2Conf/lines
```
and should be booked in 

```
Hlt/Hlt2Conf/python/Hlt2Conf/lines/bnoc/hlt2_bnoc.py
```
Note this is not necessary to run the code, it only books it with @register_line. You will notice when writing the line that the dependencies are very different. What I did was I followed the [tutorial](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/hlt2_line.html#) to write an Hlt2 line and modified it based on the run 2 line. There is a B2OC which is fully implemented with ThOr functors [here](https://gitlab.cern.ch/lhcb/Moore/-/tree/master/Hlt/Hlt2Conf/python/Hlt2Conf/lines/b_to_open_charm_thor) which is helpful as an example. 

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

# Running 
To run the code, you must write an options file. I put mine in the directory 'stack' which contains all the projects, but eventually it will go into a directory containing all the options files 
```
Hlt/Hlt2Conf/options
```
Then, follow this [tutorial](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/hlt2_line.html#running) which will tell you how to write an options file that runs over the minimum bias input. 

To run over MC files, find the file you want in bookkeeping by following [this tutorial](https://lhcb.github.io/starterkit-lessons/first-analysis-steps/bookkeeping.html). Using the path to that file, create an array in your options file which contains the paths to the files you want. The path should be 
```
root://eoslhcb.cern.ch//eos/lhcb/grid/prod/INSERT PATH YOU COPIED FROM BOOKKEEPING HERE
```
You may encounter MC files of different formats, [this tutorial](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/different_samples.html) will show you how to handle them--but it won't tell you everything! Also, some of the example options files it gives are inconsistent with what the tutorial says--for example options.input_raw_format is set in the tutorial but not in the example. I used the tutorial setting. It is also missing some things. You need to include two options not mentioned in the tutorial which are ```options.dddb_tag``` and ```options.conddb_tag```. The settings for these are specific to each decay, you can find them by following the instructions [here](https://lhcb.github.io/starterkit-lessons/first-analysis-steps/minimal-dv-job.html) under the heading "Database tags". Errors related to the decoding version can be handled following the instructions [here](https://lhcbdoc.web.cern.ch/lhcbdoc/moore/master/tutorials/different_samples.html#ft-decoding-version).

#Options File

To run over more than one file, there are a few options. At first I was running over files in eos, but I realized not all the MC files were in there. Then, I looked into using a catalog instead, follwing [these](https://lhcb.github.io/starterkit-lessons/first-analysis-steps/files-from-grid.html) instructions under the "read files remotely..." heading. An example of using inputs from eos is [her](https://github.com/umd-lhcb/hlt-umd/blob/139e95f6c6d74769452e5a81d88a84fc5feaf269/test_b_to_kpi_xdigi.py).

#Moore and DaVinci

#Moore and Ganga

I have a modified version of Moore and Rec that I need to use when submitting jobs to Ganga. For this, I use 

