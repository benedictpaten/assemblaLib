/*
 * SCAFFOLD_PATHS.h
 *
 *  Created on: 24 Feb 2011
 *      Author: benedictpaten
 */

#ifndef SCAFFOLD_PATHS_H_
#define SCAFFOLD_PATHS_H_

#include "cactus.h"
#include "sonLib.h"
#include "adjacencyClassification.h"

/*
 * Returns the length of the scaffold paths.
 */
stHash *getMaximalScaffoldPathLengths(stList *contigPaths, stList *eventStrings, CapCodeParameters *capCodeParameters);

/*
 * Gets a set of scaffold paths, each being represented as a set of contig paths.
 */
stHash *getScaffoldPaths(stList *contigPaths, stList *eventStrings, CapCodeParameters *capCodeParameters);

#endif /* SCAFFOLD_PATHS_H_ */
