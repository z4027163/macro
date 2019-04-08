#!/bin/bash

exedir=${CMSSW_BASE}/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros

melalibdir=${CMSSW_BASE}/lib/slc6_amd64_gcc630
melaincdir=${CMSSW_BASE}/src


cmsswlibdir=$CMSSW_RELEASE_BASE/lib/slc6_amd64_gcc630
cmsswincdir=$CMSSW_RELEASE_BASE/src

export LD_LIBRARY_PATH=${melalibdir}:${cmsswlibdir}:$LD_LIBRARY_PATH

if [ "$1" == "" ]; then
    echo "Please provide an arguments to the script: 4e, 4mu, 2e2mu or all"
    exit
fi

if [ $1 == "llbb_single" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_single.C HZZ4LeptonsAnalysis_llbb.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -o RunReference_llbb_single

elif [ $1 == "llbb" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C HZZ4LeptonsAnalysis_llbb.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_llbb

elif [ $1 == "btag" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C bjet_efficiency.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir}  -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_btag

elif [ $1 == "4l2b" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C HZZ4LeptonsAnalysis_4l2b.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_4l2b

elif [ $1 == "2e2mu2b" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C HZZ4LeptonsAnalysis_2e2mu2b.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_2e2mu2b

elif [ $1 == "fkj2e2mu2b" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJet_2e2mu2b.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkj2e2mu2b

elif [ $1 == "fkj2e" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJet2e.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkj2e

elif [ $1 == "fkjmu" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJet_mu.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkjmu

elif [ $1 == "fkj4mu2b" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJet_4mu.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir}  -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkj4mu2b

elif [ $1 == "fkj4e2b" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJet_4e.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkj4e2b


elif [ $1 == "fkj" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJet.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir}  -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkj

elif [ $1 == "fkj2" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJet2.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir}  -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkj2


elif [ $1 == "fkjtest" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJetTest.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkjtest

elif [ $1 == "fkjtest2" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJetTest2.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkjtest2

elif [ $1 == "fkjtestmu" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C FakeJetTest_mu.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_fkjtestmu

elif [ $1 == "rev" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C HZZ4LeptonsAnalysis_rev_4l2b.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_rev

elif [ $1 == "rev2" ]; then
    echo "Compiling $1 macros"
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C HZZ4LeptonsAnalysis_rev2_4l2b.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_rev2

elif [ $1 == "dybb" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C HZZ4LeptonsAnalysis_dybb.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_dybb
elif [ $1 == "ttbb" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C HZZ4LeptonsAnalysis_ttbb.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_ttbb
elif [ $1 == "ttbb_4m" ]; then
    echo "Compiling $1 macros" 
    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir}  -I ${cmsswincdir} compilereference_list.C HZZ4LeptonsAnalysis_ttbb_4m.C -I roccor `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir} -lCondFormatsJetMETObjects  -lJetMETCorrectionsModules -lCondFormatsBTauObjects -lCondToolsBTau  -o RunReference_ttbb_4m

fi

