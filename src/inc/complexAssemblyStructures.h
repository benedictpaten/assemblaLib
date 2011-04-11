/*
 * assemblyErrorStructures.h
 *
 *  Created on: 18 Mar 2011
 *      Author: benedictpaten
 */

#ifndef ASSEMBLYERRORSTRUCTURES_H_
#define ASSEMBLYERRORSTRUCTURES_H_

#include "cactus.h"
#include "sonLib.h"

typedef struct _CapCodeParameters {
        int32_t minimumNCount;
        int32_t maxInsertionLength;
        int32_t maxDeletionLength;
} CapCodeParameters;

CapCodeParameters *capCodeParameters_construct(int32_t minimumNCount,
                                               int32_t maxInsertionLength,
                                               int32_t maxDeletionLength);

void capCodeParameters_destruct(CapCodeParameters *capCodeParameters);

enum CapCode {
    /*
     * Switches between haplotypes
     */
    HAP_SWITCH,
    HOMO_TO_HETERO_SWITCH,
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
 * Gets a code indicating the type of structure the cap is part of.
 */
enum CapCode getCapCode(Cap *cap, stList *eventStrings, int32_t *insertLength, int32_t *deleteLength,
                        CapCodeParameters *capCodeParameters);


#endif /* ASSEMBLYERRORSTRUCTURES_H_ */
