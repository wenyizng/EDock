CC=g++
HEADER1=MathTools.h mol.h energy.h base_score.h match.h output.h translate_rotate.h MC_simulation.h REMC.h initialmol.h flex_dock.h precon.h generaterandnum.h
SOURCE1=MathTools.cpp mol.cpp energy.cpp base_score.cpp match.cpp output.cpp translate_rotate.cpp REMC1.cpp  MC_simulation2.cpp initialmol.cpp flex_dock.cpp precon.cpp generaterandnum.cpp

all: edock_flexible
edock_flexible: edock_flexible.cpp ${SOURCE1} ${HEADER1}
	${CC} ${SOURCE1} $@.cpp -o $@

clean:
	rm edock_flexible
