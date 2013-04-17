/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ASSEMBLYERRORSTRUCTURES_H_
#define ASSEMBLYERRORSTRUCTURES_H_

#include "cactus.h"
#include "sonLib.h"

/*
 * Functions to get the 'code' of an adjacency.
 */

/*
 * A wrapper structure to represent parameters to the getCapCode function.
 */
typedef struct _CapCodeParameters {
        int64_t minimumNCount;
        int64_t maxInsertionLength;
        int64_t maxDeletionLength;
} CapCodeParameters;

/*
 * Constructs a cap code parameters struct. minimumNCount is the minimum
 * number of contiguous N's to specify a scaffold gap, maxInsertionLength is the
 * maximum number size of an insertion in the chosen thread before it is reclassified as simply
 * an intra chromosomal rearrangement.
 * Similarly maxDeletionLength is the maximum size of a deletion with respect to the chosen thread before
 * it is reclassified as simply an intra chromosomal rearrangement.
 *
 */
CapCodeParameters *capCodeParameters_construct(int64_t minimumNCount,
                                               int64_t maxInsertionLength,
                                               int64_t maxDeletionLength);

/*
 * Frees the memory associated with a CapCodeParameters struct.
 */
void capCodeParameters_destruct(CapCodeParameters *capCodeParameters);

/*
 * An enum to represent the different types of adjacencies.
 * As these codes returned with respect to a given cap we call them 'capCodes'.
 */
enum CapCode {
    /*
     * Switches between haplotypes
     */
    HAP_SWITCH,
    HAP_NOTHING,

    /*
     * The end of a contig.
     */
    CONTIG_END,

    /*
     * Indications of scaffold gaps.
     */
    CONTIG_END_WITH_AMBIGUITY_GAP,
    CONTIG_END_WITH_SCAFFOLD_GAP,
    AMBIGUITY_GAP,
    SCAFFOLD_GAP,

    /*
     * The various error classes that we recognise.
     */
    ERROR_HAP_TO_HAP_SAME_CHROMOSOME,
    ERROR_HAP_TO_HAP_DIFFERENT_CHROMOSOMES,
    ERROR_HAP_TO_CONTAMINATION,
    ERROR_HAP_TO_INSERT_TO_CONTAMINATION,
    ERROR_HAP_TO_INSERT,
    ERROR_HAP_TO_INSERT_AND_DELETION,
    ERROR_HAP_TO_DELETION,
    ERROR_CONTIG_END_WITH_INSERT
};

/*
 * Gets a code indicating the type of structure the cap's adjacency is part of.
 */
enum CapCode getCapCode(Cap *cap, Cap **otherCap, stList *haplotypeEventStrings, stList *contaminationEventStrings, int64_t *insertLength, int64_t *deleteLength,
                        CapCodeParameters *capCodeParameters);



#endif /* ASSEMBLYERRORSTRUCTURES_H_ */
