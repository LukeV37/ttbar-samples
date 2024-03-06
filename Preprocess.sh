cd preprocess
root -b -q preprocess.cxx
root -q -q 'edges.cxx(0)'
root -b -q 'edges.cxx(1)'
cd ..
