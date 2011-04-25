rootPath=../../assemblaLib
cactusRootPath=${rootPath}/../cactus
include ${cactusRootPath}/include.mk

localLibPath = ../lib
libSources = impl/*.c
libHeaders = inc/*.h
cactusLibPath=${cactusRootPath}/lib

all : ${localLibPath}/assemblaLib.a

${localLibPath}/assemblaLib.a : ${libSources} ${libHeaders} ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath} -I ${cactusLibPath} -c ${libSources}
	ar rc assemblaLib.a *.o
	ranlib assemblaLib.a 
	rm *.o
	mv assemblaLib.a ${localLibPath}/
	cp ${libHeaders} ${localLibPath}/

clean : 
	rm -f *.o ${localLibPath}/*

	