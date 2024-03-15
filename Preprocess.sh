#!/bin/bash
run_ttbar=true
run_ZZ4nu=false

if [ $run_ttbar = true ] ; then
    cd preprocess
    root -b -q 'preprocess.cxx+(0)'
    root -q -q 'edges.cxx+(0)'
    root -q -q 'edges.cxx+(1)'
    mv train_data.txt ../output/ttbar
    mv train_edges.txt ../output/ttbar
    mv test_data.txt ../output/ttbar
    mv test_edges.txt ../output/ttbar
    cd ..
fi

if [ $run_ZZ4nu = true ] ; then
    cd preprocess
    root -b -q 'preprocess.cxx+(1)'
    root -q -q 'edges.cxx+(0)'
    root -q -q 'edges.cxx+(1)'
    mv train_data.txt ../output/ZZ4nu
    mv train_edges.txt ../output/ZZ4nu
    mv test_data.txt ../output/ZZ4nu
    mv test_edges.txt ../output/ZZ4nu
    cd ..
fi
