/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ASSEMBLY_STRUCTURES_H_
#define ASSEMBLY_STRUCTURES_H_

#include "cactus.h"
#include "sonLib.h"

/*
 * Basic library of functions used in tracing paths through the graph.
 */

/*
 * Gets the length of a terminal adjacency.
 */
int32_t getTerminalAdjacencyLength(Cap *cap);

/*
 * Get the sequence associated with an a terminal adjacency.
 */
char *getTerminalAdjacencySubString(Cap *cap);

/*
 * Returns non-zero iff the end contains a cap whose event is labelled with the given event string.
 */
bool hasCapInEvent(End *end, const char *eventString);

/*
 * Returns non-zero iff the end contains a cap whose event is NOT labelled with the given event string.
 */
bool hasCapNotInEvent(End *end, const char *eventString);

/*
 * Returns non-zero iff the end contains a cap with an event that labelled with one of the given event strings.
 */
bool hasCapInEvents(End *end, stList *eventStrings);

/*
 * Returns the terminal adjacency for the given cap.
 */
Cap *getTerminalCap(Cap *cap);

/*
 * Returns non-zero iff the terminal adjacency is represented in the given set of events, which are specified by the set of event strings.
 */
bool trueAdjacency(Cap *cap, stList *eventStrings);

/*
 * Returns the segment of the terminal cap.
 */
Segment *getCapsSegment(Cap *cap);

/*
 * Returns the segment of the adjacent terminal cap.
 */
Segment *getAdjacentCapsSegment(Cap *cap);

/*
 * Returns non-zero iff the caps are adjacent, initialising separation distance with the distance
 * between them.
 */
bool capsAreAdjacent(Cap *cap1, Cap *cap2, int32_t *separationDistance);

/*
 * Returns non-zero iff the ends are connected by sequences with the given events.
 */
bool endsAreConnected(End *end1, End *end2, stList *eventStrings);

/*
 * Returns non-zero iff the ends are adjacent by the sequences with the given events. If they
 * are minimumDistanceBetweenCaps is initialised to the minimum distance.
 */
bool endsAreAdjacent(End *end1, End *end2,
        int32_t *minimumDistanceBetweenCaps,
        stList *eventStrings);

#endif /* ASSEMBLY_STRUCTURES_H_ */
