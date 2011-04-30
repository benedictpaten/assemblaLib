/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef SCAFFOLD_PATHS_H_
#define SCAFFOLD_PATHS_H_

#include "cactus.h"
#include "sonLib.h"
#include "adjacencyClassification.h"

/*
 * Returns the length of the scaffold paths.
 */
stHash *getContigPathToScaffoldPathLengthsHash(stList *contigPaths, stList *eventStrings, CapCodeParameters *capCodeParameters);

/*
 * Gets a set of scaffold paths, each being represented as a set of contig paths.
 */
stHash *getScaffoldPaths(stList *contigPaths, stList *eventStrings, CapCodeParameters *capCodeParameters);

#endif /* SCAFFOLD_PATHS_H_ */
