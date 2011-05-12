/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>

#include "cactus.h"
#include "sonLib.h"
#include "contigPaths.h"
#include "adjacencyTraversal.h"
#include "adjacencyClassification.h"

int32_t getNumberOfNs(const char *string) {
    int32_t j = 0;
    for (int32_t i = 0; i < strlen(string); i++) {
        if (toupper(string[i]) == 'N') {
            j++;
        }
    }
    return j;
}

static int32_t getNumberOfNsInSegment(Segment *segment) {
    char *string = segment_getString(segment);
    int32_t i = getNumberOfNs(string);
    free(string);
    return i;
}

static int32_t getNumberOfNsInAdjacency(Cap *cap) {
    char *cA = getTerminalAdjacencySubString(cap);
    int32_t i = getNumberOfNs(cA);
    free(cA);
    return i;
}

static bool getCapGetAtEndOfPath(Cap *cap, Cap **pathEndCap,
        int32_t *pathLength, int32_t *nCount, stList *haplotypeEventStrings, stList *contaminationEventStrings) {
    //Account for length of adjacency
    *pathLength += getTerminalAdjacencyLength(cap);
    *nCount += getNumberOfNsInAdjacency(cap);

    Segment *segment = getAdjacentCapsSegment(cap);
    if (segment == NULL) {
        *pathEndCap = cap_getAdjacency(getTerminalCap(cap));
        assert(*pathEndCap != NULL);
        return 0;
    }
    Cap *adjacentCap = cap_getSide(cap) ? segment_get3Cap(segment)
    : segment_get5Cap(segment);
    assert(
            cap_getName(adjacentCap) == cap_getName(
                    cap_getAdjacency(getTerminalCap(cap))));

    End *adjacentEnd = cap_getEnd(adjacentCap);
    if (hasCapInEvents(adjacentEnd, contaminationEventStrings) || hasCapInEvents(adjacentEnd, haplotypeEventStrings)) { //hasCapNotInEvent(adjacentEnd, event_getHeader(cap_getEvent(cap)))) { //isContaminationEnd(adjacentEnd) || isHaplotypeEnd(adjacentEnd)) {
        *pathEndCap = adjacentCap;
        return 1;
    }
    *pathLength += segment_getLength(segment);
    *nCount += getNumberOfNsInSegment(segment);
    return getCapGetAtEndOfPath(cap_getOtherSegmentCap(adjacentCap),
            pathEndCap, pathLength, nCount, haplotypeEventStrings, contaminationEventStrings);
}

static int32_t getBoundingNsP(Segment *segment) {
    char *string = segment_getString(segment);
    assert(string != NULL);
    int32_t k = 0;
    for (int32_t j = 0; j < strlen(string) && j < 5; j++) {
        if (toupper(string[j]) == 'N') {
            k++;
        }
        //} else {
        //    break;
        //}
    }
    free(string);
    return k;
}

static int32_t getBoundingNs(Cap *cap) {
    assert(cap != NULL);
    Segment *segment = getCapsSegment(cap);
    if (segment == NULL) {
        return 0;
    }
    Cap *_5TerminalCap = getTerminalCap(segment_get5Cap(segment));
    Cap *_3TerminalCap = getTerminalCap(segment_get3Cap(segment));
    assert(_5TerminalCap != NULL);
    assert(_3TerminalCap != NULL);
    if (cap_getName(_5TerminalCap) == cap_getName(cap)) {
        return getBoundingNsP(segment);
    } else {
        assert(cap_getName(_3TerminalCap) == cap_getName(cap));
        return getBoundingNsP(segment_getReverse(segment));
    }
}

CapCodeParameters *capCodeParameters_construct(int32_t minimumNCount,
        int32_t maxInsertionLength, int32_t maxDeletionLength) {
    CapCodeParameters *capCodeParameters = st_malloc(sizeof(CapCodeParameters));
    capCodeParameters->minimumNCount = minimumNCount;
    capCodeParameters->maxInsertionLength = maxInsertionLength;
    capCodeParameters->maxDeletionLength = maxDeletionLength;
    return capCodeParameters;
}

void capCodeParameters_destruct(CapCodeParameters *capCodeParameters) {
    free(capCodeParameters);
}

static stSortedSet *getEventStrings(End *end, stList *eventStrings) {
    stSortedSet *eventStringsSet = stSortedSet_construct3(
            (int(*)(const void *, const void *)) strcmp, NULL);
    End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(instanceIt)) != NULL) {
        const char *header = event_getHeader(cap_getEvent(cap));
        for(int32_t i=0; i<stList_length(eventStrings); i++) {
            if(strcmp(stList_get(eventStrings, i), header) == 0) {
                stSortedSet_insert(eventStringsSet,
                                   (void *) event_getHeader(cap_getEvent(cap)));
            }
        }
    }
    end_destructInstanceIterator(instanceIt);
    return eventStringsSet;
}

static enum CapCode getHaplotypeSwitchCode(Cap *cap, stList *eventStrings) {
    Cap *adjacentCap = cap_getAdjacency(getTerminalCap(cap));
    assert(adjacentCap != NULL);
    End *end = cap_getEnd(cap);
    End *adjacentEnd = cap_getEnd(adjacentCap);
    stSortedSet *eventStringsForEnd1 = getEventStrings(end, eventStrings);
    stSortedSet *eventStringsForEnd2 = getEventStrings(adjacentEnd, eventStrings);

    assert(stSortedSet_size(eventStringsForEnd1) > 0);
    assert(stSortedSet_size(eventStringsForEnd2) > 0);

    stSortedSet *intersectionOfEventStrings = stSortedSet_getIntersection(
            eventStringsForEnd1, eventStringsForEnd2);

    enum CapCode code1 = (stSortedSet_size(intersectionOfEventStrings)
            != stSortedSet_size(eventStringsForEnd1) || stSortedSet_size(
                    intersectionOfEventStrings)
            != stSortedSet_size(eventStringsForEnd2)) ? HAP_SWITCH
    : HAP_NOTHING;

    stSortedSet_destruct(eventStringsForEnd1);
    stSortedSet_destruct(eventStringsForEnd2);
    stSortedSet_destruct(intersectionOfEventStrings);

    return code1;
}

enum CapCode getCapCode(Cap *cap, stList *haplotypeEventStrings, stList *contaminationEventStrings, int32_t *insertLength,
        int32_t *deleteLength, CapCodeParameters *capCodeParameters) {
    assert(hasCapInEvents(cap_getEnd(cap), haplotypeEventStrings));
    if (trueAdjacency(cap, haplotypeEventStrings)) {
        return getHaplotypeSwitchCode(cap, haplotypeEventStrings);
    }
    *insertLength = 0;
    *deleteLength = 0;

    End *end = cap_getEnd(cap);
    Cap *pathEndCap = NULL;
    int32_t pathLength = 0, nCount = 0;
    bool pathEndsOnStub = !getCapGetAtEndOfPath(cap, &pathEndCap, &pathLength,
            &nCount, haplotypeEventStrings, contaminationEventStrings);
    assert(pathLength >= 0);
    assert(nCount >= 0);
    assert(pathEndCap != NULL);
    End *otherPathEnd = cap_getEnd(pathEndCap);
    nCount += getBoundingNs(cap) + getBoundingNs(pathEndCap);

    if (pathEndsOnStub) {
        assert(!hasCapInEvents(otherPathEnd, contaminationEventStrings)); //Can not test for hap event strings, as a stub end may contain the reference.
        //assert(!hasCapInEvents(otherPathEnd, haplotypeEventStrings)); //isHaplotypeEnd(otherPathEnd) && !isContaminationEnd(otherPathEnd));
        return pathLength == 0 ? CONTIG_END
        : (nCount >= 1 ? (nCount >= capCodeParameters->minimumNCount ? CONTIG_END_WITH_SCAFFOLD_GAP
                        : CONTIG_END_WITH_AMBIGUITY_GAP)
                : ERROR_CONTIG_END_WITH_INSERT);
    }

    if (hasCapInEvents(otherPathEnd, haplotypeEventStrings)) {
        if (endsAreConnected(end, otherPathEnd, haplotypeEventStrings)) {
            int32_t minimumHaplotypeDistanceBetweenEnds;
            //Establish if indel or order breaking rearrangement
            if (endsAreAdjacent(end, otherPathEnd,
                            &minimumHaplotypeDistanceBetweenEnds, haplotypeEventStrings)) {
                *insertLength = pathLength;
                *deleteLength = minimumHaplotypeDistanceBetweenEnds;
                if (nCount >= capCodeParameters->minimumNCount) { //Insertion was scaffold gap
                    return SCAFFOLD_GAP;
                }
                if (nCount >= 1) {
                    return AMBIGUITY_GAP;
                }
                if (pathLength > 0) { //Had insertion
                    if (minimumHaplotypeDistanceBetweenEnds > 0) { //Had deletion also
                        return (minimumHaplotypeDistanceBetweenEnds
                                >= capCodeParameters->maxDeletionLength
                                || pathLength
                                >= capCodeParameters->maxInsertionLength) ? ERROR_HAP_TO_HAP_SAME_CHROMOSOME
                        : ERROR_HAP_TO_INSERT_AND_DELETION;
                    }
                    return pathLength >= capCodeParameters->maxInsertionLength ? ERROR_HAP_TO_HAP_SAME_CHROMOSOME
                    : ERROR_HAP_TO_INSERT;
                } else { //Deletion, possibly with Ns
                    assert(minimumHaplotypeDistanceBetweenEnds > 0);
                    return minimumHaplotypeDistanceBetweenEnds
                    >= capCodeParameters->maxDeletionLength ? ERROR_HAP_TO_HAP_SAME_CHROMOSOME
                    : ERROR_HAP_TO_DELETION;
                }
            } else {
                return ERROR_HAP_TO_HAP_SAME_CHROMOSOME;
            }
        } else {
            return ERROR_HAP_TO_HAP_DIFFERENT_CHROMOSOMES;
        }
    } else {
        return pathLength == 0 ? ERROR_HAP_TO_CONTAMINATION
        : ERROR_HAP_TO_INSERT_TO_CONTAMINATION;
    }
}
