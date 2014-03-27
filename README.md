BtagCode
=======================

1) First thing to do is:

    cd CMSSW_5_X_Y/src
    cmsenv
    
note that this package works only in a CMSSW_5_X_Y

2) In BtagCode directory:

    source scripts/setup.sh
    
3) Compile:

    make 
    make exe
    
4) Run:

    ./bin/Make_BtagEfficiencyMaps.exe cfg/Make_BtagEfficiencyMaps.cfg
    
5) NOTE: 
   in the cfg file there are:

    inputList = list of file to run
    inputTree = name of the tree (all the files must have the same tree name)
    outputName = output name
    

BtagCode on lxbatch
==================================

1) First thing to do is:

    cd CMSSW_5_X_Y/src
    cmsenv
    
note that this package works only in a CMSSW_5_X_Y

2) In BtagCode directory, compile:

    make 
    make exe
    
3) In BtagCode/test/lxbatch directory prepare the jobs:

    perl launchJobs_lxbatch_eos.pl params_lxbatch.CFG
    
    params_lxbatch.CFG has the following input parameters:
    
    - BASEDir: complete path of this lxbatch directory, eg:
    
        /afs/cern.ch/work/b/bmarzocc/BtagCode/test/lxbatch
        
    - CMSSWPath: CMSSW_5_X_Y release path from which set up the environment, eg:
        
        /afs/cern.ch/work/b/bmarzocc/CMSSW_5_3_10/src/
        
    - ProgramName: executable file name, eg:
    
        Make_BtagEfficiencyMaps
        
    - JOBCfgTemplate: cfg file to run, USE THE TEMPLATE:
    
        /afs/cern.ch/work/b/bmarzocc/BtagCode/test/lxbatch/parserParams_template.cfg
        
    - LISTOFSamples: txt file of the list of directories that contain the root files, eg of path into the txt:
      
        /store/cmst3/user/obondu/H2GGLOBE/Radion/reduced/radion_reduction_v11/mc RadionToHH_2Gamma_2b_M-1100_TuneZ2star_8TeV-Madgraph_pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1

        where the directory path and the directory have to be separated by a spacetab
        
    - OUTPUTSAVEPath: directory where to save the output files (also a eos directory), eg:
         
        /afs/cern.ch/work/b/bmarzocc/BtagCode/test/lxbatch/
        
    - OUTPUTFILEName: name of a single job output root file, eg:
        
        btagEfficiencies_Radion1100
    
    - JOBModulo: numeber of input file run per job
    
    - Queue
    
3) Launch the jobs (file lancia.sh automatically created):

    sh lancia.sh

4) Add all the final maps in a unique file:

    source scripts/setup.sh
    ./bin/AddHistos2D.exe cfg/AddHistos.cfg
    
in the cfg file there are:
    
    inputDir = the directory where the .root files to be added are stored
    outputName = name of the output file
    



