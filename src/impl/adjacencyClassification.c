/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>

#include "cactus.h"
#include "sonLib.h"
#include "haplotypePaths.h"
#include "adjacencyTraversal.h"
#include "adjacencyClassification.h"

static bool getCapGetAtEndOfPath(Cap *cap, Cap **pathEndCap, stList *segments) {
    Segment *segment = getAdjacentCapsSegment(cap);
    if (segment == NULL) {
        *pathEndCap = cap_getAdjacency(getTerminalCap(cap));
        assert(*pathEndCap != NULL);
        return 0;
    }
    Cap *adjacentCap = cap_getSide(cap) ? segment_get3Cap(segment) : segment_get5Cap(segment);
    assert(cap_getName(adjacentCap) == cap_getName(cap_getAdjacency(
                            getTerminalCap(cap))));
    End *adjacentEnd = cap_getEnd(adjacentCap);
    if (hasCapNotInEvent(adjacentEnd, event_getHeader(cap_getEvent(cap)))) { //isContaminationEnd(adjacentEnd) || isHaplotypeEnd(adjacentEnd)) {
        *pathEndCap = adjacentCap;
        return 1;
    }
    if (segments != NULL) {
        stList_append(segments, segment);
    }
    return getCapGetAtEndOfPath(cap_getOtherSegmentCap(adjacentCap), pathEndCap, segments);
}

static int32_t getNumberOfNs(Segment *segment) {
    char *string = segment_getString(segment);
    int32_t j = 0;
    for (int32_t i = 0; i < strlen(string); i++) {
        if (toupper(string[i]) == 'N') {
            j++;
        }
    }
    free(string);
    return j;
}

static int32_t getBoundingNsP(Segment *segment) {
    char *string = segment_getString(segment);
    assert(string != NULL);
    int32_t k = 0;
    for (int32_t j = 0; j < strlen(string); j++) {
        if (string[j] == 'N' || string[j] == 'n') {
            k++;
        } else {
            break;
        }
    }
    free(string);
    return k;
}

static int32_t getBoundingNs(Cap *cap) {
    assert(cap != NULL);
    Segment *segment = getAdjacentCapsSegment(cap);
    if (segment == NULL) {
        return 0;
    }
    Cap *_5TerminalCap = getTerminalCap(segment_get5Cap(segment));
    Cap *_3TerminalCap = getTerminalCap(segment_get3Cap(segment));
    assert(_5TerminalCap != NULL);
    assert(_3TerminalCap != NULL);
    if (cap_getName(cap_getAdjacency(_5TerminalCap)) == cap_getName(cap)) {
        return getBoundingNsP(segment);
    } else {
        assert(cap_getName(cap_getAdjacency(_3TerminalCap)) == cap_getName(cap));
        return getBoundingNsP(segment_getReverse(segment));
    }
}

CapCodeParameters *capCodeParameters_construct(int32_t minimumNCount,
                                               int32_t maxInsertionLength,
                                               int32_t maxDeletionLength) {
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
    stSortedSet *eventStringsSet = stSortedSet_construct3((int (*)(const void *, const void *))strcmp, NULL);
    End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
    Cap *cap;
    while((cap = end_getNext(instanceIt)) != NULL) {
        stSortedSet_insert(eventStringsSet, (void *)event_getHeader(cap_getEvent(cap)));
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

    stSortedSet *intersectionOfEventStrings = stSortedSet_getIntersection(eventStringsForEnd1, eventStringsForEnd2);

    enum CapCode code1 = (stSortedSet_size(intersectionOfEventStrings) != stSortedSet_size(eventStringsForEnd1) || stSortedSet_size(intersectionOfEventStrings) != stSortedSet_size(eventStringsForEnd2)) ? HAP_SWITCH : HAP_NOTHING;

    stSortedSet_destruct(eventStringsForEnd1);
    stSortedSet_destruct(eventStringsForEnd2);
    stSortedSet_destruct(intersectionOfEventStrings);

    return code1;
}

enum CapCode getCapCode(Cap *cap, stList *eventStrings, int32_t *insertLength, int32_t *deleteLength,
        CapCodeParameters *capCodeParameters) {
    assert(hasCapInEvents(cap_getEnd(cap), eventStrings));
    if (trueAdjacency(cap, eventStrings)) {
    	return getHaplotypeSwitchCode(cap, eventStrings);
    }
    *insertLength = 0;
    *deleteLength = 0;
    stList *segments = stList_construct();

    End *end = cap_getEnd(cap);
    Cap *pathEndCap = NULL;
    bool pathEndsOnStub = !getCapGetAtEndOfPath(cap, &pathEndCap, segments);
    assert(pathEndCap != NULL);
    End *otherPathEnd = cap_getEnd(pathEndCap);

    bool hadSegments = stList_length(segments) > 0;
    int32_t nCount = getBoundingNs(cap) + getBoundingNs(pathEndCap);
    int32_t totalSegmentLength = 1;
    for (int32_t i = 0; i < stList_length(segments); i++) {
        nCount += getNumberOfNs(stList_get(segments, i));
        totalSegmentLength += segment_getLength(stList_get(segments, i));
    }
    stList_destruct(segments);
    if (pathEndsOnStub) {
        assert(!hasCapInEvents(otherPathEnd, eventStrings)); //isHaplotypeEnd(otherPathEnd) && !isContaminationEnd(otherPathEnd));
        return !hadSegments ? CONTIG_END : (nCount >= 1 ? (nCount >= capCodeParameters->minimumNCount ? CONTIG_END_WITH_SCAFFOLD_GAP
                : CONTIG_END_WITH_AMBIGUITY_GAP) : ERROR_CONTIG_END_WITH_INSERT);
    }

    //assert(hasCapInEvents(cap_getEnd(cap), events)); //isHaplotypeEnd(otherPathEnd) || isContaminationEnd(otherPathEnd));
    //const char *hapNames[2] = { "hapA1", "hapA2" };
    if (hasCapInEvents(otherPathEnd, eventStrings)) { //isHaplotypeEnd(otherPathEnd)) {
        if (endsAreConnected(end, otherPathEnd, eventStrings)) {
            int32_t minimumHaplotypeDistanceBetweenEnds;
            //Establish if indel or order breaking rearrangement
            if (endsAreAdjacent(end, otherPathEnd, &minimumHaplotypeDistanceBetweenEnds, eventStrings)) {
                *insertLength = totalSegmentLength;
                *deleteLength = minimumHaplotypeDistanceBetweenEnds;
                if (nCount >= capCodeParameters->minimumNCount) { //Insertion was scaffold gap
                    return SCAFFOLD_GAP;
                }
                if (nCount >= 1) {
                    return AMBIGUITY_GAP;
                }
                if (hadSegments) { //Had insertion
                    if (minimumHaplotypeDistanceBetweenEnds > 0) { //Had deletion also
                        return (minimumHaplotypeDistanceBetweenEnds >= capCodeParameters->maxDeletionLength || totalSegmentLength >= capCodeParameters->maxInsertionLength) ? ERROR_HAP_TO_HAP_SAME_CHROMOSOME : ERROR_HAP_TO_INSERT_AND_DELETION;
                    }
                    return totalSegmentLength >= capCodeParameters->maxInsertionLength ? ERROR_HAP_TO_HAP_SAME_CHROMOSOME : ERROR_HAP_TO_INSERT;
                } else { //Deletion, possibly with Ns
                    assert(minimumHaplotypeDistanceBetweenEnds > 0);
                    return minimumHaplotypeDistanceBetweenEnds >= capCodeParameters->maxDeletionLength ? ERROR_HAP_TO_HAP_SAME_CHROMOSOME : ERROR_HAP_TO_DELETION;
                }
            } else {
                return ERROR_HAP_TO_HAP_SAME_CHROMOSOME;
            }
        } else {
            return ERROR_HAP_TO_HAP_DIFFERENT_CHROMOSOMES;
        }
    } else {
        return !hadSegments ? ERROR_HAP_TO_CONTAMINATION : ERROR_HAP_TO_INSERT_TO_CONTAMINATION;
    }
}
