j = Job( name = 'MC_job')
myApp = GaudiExec()
myApp.directory='/afs/cern.ch/work/e/ejiang/public/Moore_nightly'
j.application = myApp
j.application.platform="x86_64_v2-centos7-gcc11-opt"
j.application.options = [
    "/afs/cern.ch/work/e/ejiang/public/Moore_nightly/test_b_to_kpi_xdigi.py"]

j.application.readInputData("/afs/cern.ch/work/e/ejiang/public/stack2/MC_Upgrade_12101401_Beam7000GeVUpgradeMagDownNu7.625nsPythia8_Sim10Up08_XDIGI.py")

j.backend = Dirac()
j.outputfiles = [ LocalFile('test_b_to_kpi_MC_with_tck.ldst'), LocalFile('my_test_b_to_kpi_MC.tck.json') ]
j.splitter = SplitByFiles(filesPerJob=100)
j.backend.settings['CPUTime'] = 100000
j.submit()
