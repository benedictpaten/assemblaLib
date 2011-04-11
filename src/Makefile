rootPath=../../git/cactus/
include ${rootPath}include.mk

libSources = impl/*.c
libHeaders = inc/*.h 

all : ${libPath}/alignmentPathLib.a

${libPath}/alignmentPathLib.a : ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc alignmentPathLib.a *.o
	ranlib alignmentPathLib.a 
	rm *.o
	mv alignmentPathLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/

clean : 
	rm -f *.o ${libPath}/alignmentPathLib.a

	