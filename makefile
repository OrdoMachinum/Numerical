ifdef OS
	OUT = bin/numplay.exe
else
	OUT = bin/numplay
endif

$(OUT) :
	g++ -o $(OUT) numPlay.cpp GuassNewton.cpp csvtools.cpp