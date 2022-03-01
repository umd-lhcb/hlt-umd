j = Job( name = 'test_job')
myApp = GaudiExec()
myApp.directory='/afs/cern.ch/work/e/ejiang/public/Moore_nightly'
j.application = myApp

j.application.platform="x86_64_v2-centos7-gcc11-opt"
j.application.options = [
    "/afs/cern.ch/work/e/ejiang/public/Moore_nightly/test_b_to_kpi.py"]
j.application.readInputData("/afs/cern.ch/work/e/ejiang/public/stack3/MC_Upgrade_Beam7000GeVUpgradeMagDownNu7.625nsPythia8_Sim09cUp02_RecoUp01_Trig0x52000000_30000000_LDST.py")
j.backend = Dirac()
j.outputfiles = [ LocalFile('test_b_to_kpi_minbias_with_tck.LDST'), LocalFile('my_test_b_to_kpi_minbias.tck.json') ]
j.splitter = SplitByFiles(filesPerJob=100)
j.backend.settings['CPUTime'] = 100000
j.submit()
