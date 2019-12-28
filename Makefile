CC=g++
FLAGS=-std=c++11
LINK=-lgsl -lgslcblas
MAGMA_FLAGS := -DADD_ -I/usr/local/magma/include -I/usr/local/cuda-10.0/include
MAGMA_LIBS := -L/usr/local/magma/lib -L/usr/local/cuda-10.0/lib64 -L/usr/lib64/openmpi/lib/ -fopenmp -lmagma  -lcublas -lcudart -lcusparse
LAPACK_LIBS := -L/usr/lib64 -L/usr/lib64/atlas -llapack -lpthread -ltatlas

all: pmaker efs dos netrlax psi6 lrat writepoly smap dos_avg adisp prat\
	proj_calc svg_write mkbin nencode extract stlwrite

pmaker:
	$(CC) $(FLAGS) src/drawpoints.cpp -o pmaker $(LINK) -pthread

efs:
	$(CC) $(FLAGS) $(MAGMA_FLAGS) -o efs src/eig_freqs.cpp $(LAPACK_LIBS)\
        $(MAGMA_LIBS)

dos:
	$(CC) $(FLAGS) -o dos src/freq_dos.cpp

netrelax:
	$(CC) $(FLAGS) -o netrelax src/relax_network.cpp $(LINK) -pthread

psi6:
	$(CC) $(FLAGS) -o psi6 src/psi_6.cpp $(LINK)

lrat:
	$(CC) $(FLAGS) -lm -o lrat src/length_ratio_calc.cpp

writepoly:
	$(CC) $(FLAGS) -lm -o writepoly src/write_poly_file.cpp

smap:
	$(CC) $(FLAGS) -lm -o smap src/make_strain_map.cpp

dos_avg:
	$(CC) $(FLAGS) -o dos_avg src/dos_average.cpp

adisp:
	$(CC) $(FLAGS) -o adisp src/affine_displace.cpp $(LINK)

prat:
	$(CC) $(FLAGS) -o prat src/calc_prat.cpp

proj_calc:
	$(CC) $(FLAGS) -o proj_calc src/projection_calc.cpp

svg_write:
	$(CC) $(FLAGS) -o svg_write src/write_svg.cpp
mkbin:
	$(CC) $(FLAGS) -o mkbin src/make_binary.cpp
nencode:
	$(CC) $(FLAGS) -o nencode src/encode_network.cpp

extract:
	$(CC) $(FLAGS) -o extract src/extract_range.cpp

src/triangle.o:
	gcc -O -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib -DTRILIBRARY -c -o src/triangle.o src/triangle.c

stlwrite: src/triangle.o
	$(CC) $(FLAGS) -o stlwrite src/write_stl.cpp src/triangle.o -lm

trihier: src/triangle.o
	$(CC) $(FLAGS) -o trihier src/triangulated_hierarchy.cpp src/VoronoiDiagramGenerator.cpp src/triangle.o -lm $(LINK)

nwh:
	$(CC) $(FLAGS) -o nwh src/network_with_hole.cpp -lm

mdm:
	$(CC) $(FLAGS) -o mdm src/make_deviation_map.cpp

cla:
	$(CC) $(FLAGS) -o cla src/cut_line_analysis.cpp 

pmclean:
	rm pmaker

eigclean:
	rm efs

dosclean:
	rm dos

relclean:
	rm netrelax

psi6clean:
	rm psi6

lratclean:
	rm lrat

polyclean:
	rm writepoly

smclean:
	rm smap

davg_clean:
	rm dos_avg

adclean:
	rm adisp

ratclean:
	rm prat

projclean:
	rm proj_calc

svgclean:
	rm svg_write
binclean:
	rm mkbin
enclean:
	rm nencode
exclean:
	rm extract

stlclean:
	rm stlwrite

triclean:
	rm src/triangle.o

thclean:
	rm trihier

hclean:
	rm nwh

dclean:
	rm mdm

clclean:
	rm cla

purge:
	rm pmaker
	rm efs
	rm dos
	rm netrlax
	rm psi6
	rm lrat
	rm writepoly
	rm smap
	rm dos_avg 
	rm adisp 
	rm prat
	rm proj_calc
	rm svg_write
	rm mkbin
	rm nencode
	rm extract
	rm stlwrite
	rm trihier
	rm nwh
	rm mdm
	rm cla
