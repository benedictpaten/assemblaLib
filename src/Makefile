rootPath=../../cactus/
include ${rootPath}include.mk

libSources = impl/*.c
libHeaders = inc/*.h 

all : ${libPath}/assemblaLib.a

${libPath}/assemblaLib.a : ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc assemblaLib.a *.o
	ranlib assemblaLib.a 
	rm *.o
	mv assemblaLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/

clean : 
	rm -f *.o ${libPath}/assemblaLib.a

	