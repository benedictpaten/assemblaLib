rootPath=../../cactus/
include ${rootPath}include.mk

localLibPath = ../lib
libSources = impl/*.c
libHeaders = inc/*.h 

all : ${localLibPath}/assemblaLib.a

${localLibPath}/assemblaLib.a : ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc assemblaLib.a *.o
	ranlib assemblaLib.a 
	rm *.o
	mv assemblaLib.a ${localLibPath}/
	cp ${libHeaders} ${localLibPath}/

clean : 
	rm -f *.o ${localLibPath}/*

	