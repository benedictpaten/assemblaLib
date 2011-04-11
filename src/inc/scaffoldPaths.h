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

/*
 * Returns the length of the scaffold paths.
 */
stHash *getMaximalScaffoldPathLengths(stList *haplotypePaths, stList *eventStrings, CapCodeParameters *capCodeParameters);

/*
 * Gets a set of scaffold paths, each being represented as a set of haplotype paths.
 */
stHash *getScaffoldPaths(stList *haplotypePaths, stList *eventStrings, CapCodeParameters *capCodeParameters);

#endif /* SCAFFOLD_PATHS_H_ */
