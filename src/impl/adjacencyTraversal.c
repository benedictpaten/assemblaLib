/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "cactus.h"
#include "adjacencyTraversal.h"

bool getTerminalAdjacencyLength_ignoreAdjacencies = 0;

char *getTerminalAdjacencySubString(Cap *cap) {
    if(getTerminalAdjacencyLength_ignoreAdjacencies) {
        return stString_copy("");
    }
    cap = getTerminalCap(cap);
    cap = cap_getStrand(cap) ? cap : cap_getReverse(cap); //This ensures the asserts are as expected.
    Cap *adjacentCap = cap_getAdjacency(cap);
    int32_t i = cap_getCoordinate(cap) - cap_getCoordinate(adjacentCap);
    assert(i != 0);
    if (i > 0) {
        assert(cap_getSide(cap));
        assert(!cap_getSide(adjacentCap));
        return sequence_getString(cap_getSequence(cap),
                cap_getCoordinate(adjacentCap) + 1, i - 1, 1);
    } else {
        assert(cap_getSide(adjacentCap));
        assert(!cap_getSide(cap));
        return sequence_getString(cap_getSequence(cap), cap_getCoordinate(cap) + 1, -i - 1, 1);
    }
}

int32_t getTerminalAdjacencyLength(Cap *cap) {
    if(getTerminalAdjacencyLength_ignoreAdjacencies) {
        return 0;
    }
    char *cA = getTerminalAdjacencySubString(cap);
    int32_t i = strlen(cA);
    free(cA);
    return i;
}

Cap *getTerminalCap(Cap *cap) {
    Flower *nestedFlower = group_getNestedFlower(end_getGroup(cap_getEnd(cap)));
    if (nestedFlower != NULL) {
        Cap *nestedCap = flower_getCap(nestedFlower, cap_getName(cap));
        assert(nestedCap != NULL);
        return getTerminalCap(cap_getOrientation(cap) ? nestedCap : cap_getReverse(nestedCap));
    }
    return cap;
}

Segment *getCapsSegment(Cap *cap) {
    if (cap_getSegment(cap) != NULL) {
        return cap_getSegment(cap);
    }
    assert(!end_isBlockEnd(cap_getEnd(cap)));
    assert(end_isStubEnd(cap_getEnd(cap)));
    //Walk up to get the next adjacency.
    Group *parentGroup = flower_getParentGroup(end_getFlower(cap_getEnd(cap)));
    if (parentGroup != NULL) {
        Cap *parentCap = flower_getCap(group_getFlower(parentGroup), cap_getName(cap));
        if (parentCap != NULL) {
            assert(cap_getOrientation(parentCap));
            if (!cap_getOrientation(cap)) {
                parentCap = cap_getReverse(parentCap);
            }
            return getCapsSegment(parentCap);
        } else { //Cap must be a free stub end.
            assert(0); //Not in the current alignments.
            assert(end_isFree(cap_getEnd(cap)));
        }
    }
    return NULL;
}

Segment *getAdjacentCapsSegment(Cap *cap) {
    cap = getTerminalCap(cap);
    cap = cap_getAdjacency(cap);
    assert(cap != NULL);
    return getCapsSegment(cap);
}

static bool capHasGivenEvents(Cap *cap, stList *eventStrings) {
    const char *headerSeq = event_getHeader(cap_getEvent(cap));
    for (int32_t i = 0; i < stList_length(eventStrings); i++) {
        if (strcmp(headerSeq, stList_get(eventStrings, i)) == 0) {
            return 1;
        }
    }
    return 0;
}

bool trueAdjacency(Cap *cap, stList *eventStrings) {
	if(getTerminalAdjacencyLength(cap) > 0) {
        return 0;
    }
	cap = getTerminalCap(cap);
    assert(cap != NULL);
    Cap *otherCap = cap_getAdjacency(cap);
    assert(otherCap != NULL);
    assert(cap_getAdjacency(otherCap) == cap);
    //So is the adjacency present in one of the haplotypes? That's what we're going to answer..
    End *otherEnd = end_getPositiveOrientation(cap_getEnd(otherCap));
    End_InstanceIterator *endInstanceIt = end_getInstanceIterator(cap_getEnd(cap));
    Cap *cap2;
    while ((cap2 = end_getNext(endInstanceIt)) != NULL) {
        Cap *otherCap2 = cap_getAdjacency(cap2);
        assert(otherCap2 != NULL);
        if (otherEnd == end_getPositiveOrientation(cap_getEnd(otherCap2))) {
            //const char *eventName = event_getHeader(cap_getEvent(cap2));
            assert(event_getHeader(cap_getEvent(cap2)) == event_getHeader(
                            cap_getEvent(otherCap2)));
            if (capHasGivenEvents(cap2, eventStrings)) { //strcmp(eventName, "hapA1") == 0 || strcmp(eventName, "hapA2") == 0) {
                if(getTerminalAdjacencyLength(cap2) == 0) {
                    end_destructInstanceIterator(endInstanceIt);
                    return 1;
                }
            }
        }
    }
    end_destructInstanceIterator(endInstanceIt);
    return 0;
}

bool hasCapInEvent(End *end, const char *eventString) {
    Cap *cap;
    End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
    while ((cap = end_getNext(instanceIt)) != NULL) {
        if (strcmp(event_getHeader(cap_getEvent(cap)), eventString) == 0) {
            end_destructInstanceIterator(instanceIt);
            return 1;
        }
    }
    end_destructInstanceIterator(instanceIt);
    return 0;
}

bool hasCapNotInEvent(End *end, const char *eventString) {
    Cap *cap;
    End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
    while ((cap = end_getNext(instanceIt)) != NULL) {
        if (strcmp(event_getHeader(cap_getEvent(cap)), eventString) != 0) {
            end_destructInstanceIterator(instanceIt);
            return 1;
        }
    }
    end_destructInstanceIterator(instanceIt);
    return 0;
}

bool hasCapInEvents(End *end, stList *eventStrings) {
    for(int32_t i=0; i<stList_length(eventStrings); i++) {
        if(hasCapInEvent(end, stList_get(eventStrings, i))) {
            return 1;
        }
    }
    return 0;
}

bool endsAreConnected(End *end1, End *end2, stList *eventStrings) {
    if (end_getName(end1) == end_getName(end2)) { //Then the ends are the same and are part of the same chromosome by definition.
        End_InstanceIterator *instanceIterator = end_getInstanceIterator(end1);
        Cap *cap1;
        while ((cap1 = end_getNext(instanceIterator)) != NULL) {
            if (capHasGivenEvents(cap1, eventStrings)) {
                end_destructInstanceIterator(instanceIterator);
                return 1;
            }
        }
        return 0;
    }
    End_InstanceIterator *instanceIterator = end_getInstanceIterator(end1);
    Cap *cap1;
    while ((cap1 = end_getNext(instanceIterator)) != NULL) {
        if (capHasGivenEvents(cap1, eventStrings)) {
            End_InstanceIterator *instanceIterator2 = end_getInstanceIterator(end2);
            Cap *cap2;
            while ((cap2 = end_getNext(instanceIterator2)) != NULL) {
                assert(cap_getName(cap2) != cap_getName(cap1)); //This could only happen if end1 == end2
                if (sequence_getMetaSequence(cap_getSequence(cap1)) == sequence_getMetaSequence(cap_getSequence(cap2))) {
                    assert(strcmp(event_getHeader(cap_getEvent(cap1)),
                                    event_getHeader(cap_getEvent(cap2))) == 0);
                    assert(cap_getPositiveOrientation(cap1)
                            != cap_getPositiveOrientation(cap2));
                    assert(cap_getName(cap1) != cap_getName(cap2));
                    //they could have the same coordinate if they represent two ends of a block of length 1.

                    end_destructInstanceIterator(instanceIterator);
                    end_destructInstanceIterator(instanceIterator2);
                    return 1;
                }
            }
            end_destructInstanceIterator(instanceIterator2);
        }
    }
    end_destructInstanceIterator(instanceIterator);
    return 0;
}

bool capsAreAdjacent(Cap *cap1, Cap *cap2, int32_t *separationDistance) {
    if (cap_getName(cap2) != cap_getName(cap1) && cap_getCoordinate(cap1) != cap_getCoordinate(cap2)) { //This can happen if end1 == end2
        if (sequence_getMetaSequence(cap_getSequence(cap1)) == sequence_getMetaSequence(cap_getSequence(cap2))) {
            assert(strcmp(event_getHeader(cap_getEvent(cap1)), event_getHeader(
                                    cap_getEvent(cap2))) == 0);
            assert(cap_getPositiveOrientation(cap1)
                    != cap_getPositiveOrientation(cap2));
            assert(cap_getName(cap1) != cap_getName(cap2));
            assert(sequence_getMetaSequence(cap_getSequence(cap1))
                                == sequence_getMetaSequence(cap_getSequence(cap2)));

            if (!cap_getStrand(cap1)) {
                cap1 = cap_getReverse(cap1);
            }
            if (!cap_getStrand(cap2)) {
                cap2 = cap_getReverse(cap2);
            }
            assert(cap_getStrand(cap1));
            assert(cap_getStrand(cap2));
            if (cap_getCoordinate(cap1) < cap_getCoordinate(cap2)) {
                if (!cap_getSide(cap1) && cap_getSide(cap2)) {
                    *separationDistance = cap_getCoordinate(cap2) - cap_getCoordinate(cap1) - 1; //The minus 1, to give the length of the sequence between the two caps.
                    return 1;
                }
            } else {
                if (cap_getSide(cap1) && !cap_getSide(cap2)) {
                    *separationDistance = cap_getCoordinate(cap1) - cap_getCoordinate(cap2) - 1;
                    return 1;
                }
            }
        }
    }
    return 0;
}

bool endsAreAdjacent2(End *end1, End *end2, Cap **returnCap1, Cap **returnCap2, int32_t *minimumDistanceBetweenHaplotypeCaps, stList *eventStrings) {
    End_InstanceIterator *instanceIterator = end_getInstanceIterator(end1);
    Cap *cap1;
    *returnCap1 = NULL;
    *returnCap2 = NULL;
    *minimumDistanceBetweenHaplotypeCaps = INT32_MAX;
    bool areAdjacent = 0;
    while ((cap1 = end_getNext(instanceIterator)) != NULL) {
        if (capHasGivenEvents(cap1, eventStrings)) {
            End_InstanceIterator *instanceIterator2 = end_getInstanceIterator(end2);
            Cap *cap2;
            while ((cap2 = end_getNext(instanceIterator2)) != NULL) {
                int32_t i;
                if (capsAreAdjacent(cap1, cap2, &i)) {
                    areAdjacent = 1;
                    if (i < *minimumDistanceBetweenHaplotypeCaps) {
                        *minimumDistanceBetweenHaplotypeCaps = i;
                        *returnCap1 = cap1;
                        *returnCap2 = cap2;
                    }
                }
            }
            end_destructInstanceIterator(instanceIterator2);
        }
    }
    end_destructInstanceIterator(instanceIterator);
    return areAdjacent;
}

bool endsAreAdjacent(End *end1, End *end2, int32_t *minimumDistanceBetweenHaplotypeCaps, stList *eventStrings) {
    Cap *cap1 = NULL, *cap2 = NULL;
    return endsAreAdjacent2(end1, end2, &cap1, &cap2, minimumDistanceBetweenHaplotypeCaps, eventStrings);
}
