rm lib
rm include
ln -s /home/issp/materiapps/tool/Komega/Komega-2.0.0-0/lib/ .
ln -s /home/issp/materiapps/tool/Komega/Komega-2.0.0-0/include/ .
source /home/issp/materiapps/tool/Komega/Komegavars.sh 

icc -g Komega.c  -DDSFMT_MEXP=19937 -DHAVE_SSE2 dSFMT.c -qopenmp -O3 -o Komega.out -L./lib -lkomega -I./include  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread 
