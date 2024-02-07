run_ttbar=true
run_ZZ4nu=true

cd ./refine
make

if [ $run_ttbar = true ] ; then
    echo "Refining ttbar..."
    ./d_ana ../data/user.khanov.mc15_14TeV.600012.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.r12573_mc_trk_Akt4EMTo/*
    mv refined.root refined_ttbar.root
    echo "Done Refining ttbar..."
fi

if [ $run_ZZ4nu = true ] ; then
    echo "Refining ZZ4nu..."
    ./d_ana ../data/user.khanov.mc15.mc15_14TeV.600026.PhH7EG_NNPDF3_AZNLO_VBFH125_ZZ4nu_MET75.r13618_mc_trk_Akt4EMTo/*
    mv refined.root refined_ZZ4nu.root
    echo "Done Refining ZZ4nu..."
fi

#make clean
cd ..
echo "DONE :)"
