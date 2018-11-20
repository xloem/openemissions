#pragma once

#include "types.hpp"

struct RecBufMeta
{
  Scalar freq;
  Scalar rate;
  uint64_t sampleTime;
};

/*
 * - [ ] I'm leaving samples out of the accumulator.  that's inaccurate, a little, but likely ok
 *    -> is it inaccurate only because of sample loss? or is there more?
 * - [ ] I could include these samples in the accumulator.  I could add them manually, sort the matrix by
 *   row length, or have a different matrix for every period length.
 * - [X] wrt _shortest, will this work with merging many accumulated waveforms together?
 *     if some of them are short, will this cause data loss in the longer ones?
 *   -> habit is to reconsider _shortest at every incoming buffer
 *      so when adding, it should adopt value of new waveform to match habit
 *
 *      really we're compensting for inaccuracies here.  the assumption is that later periods are
 *      more like what the future will hold, and have been more accurately decided upon.
 *      good enough for now
 *
 * I'd like to add _all_ the samples in the accumulator, but since I'll need a new kind of algorithm for that,
 * maybe I'll see if it's necessary first.
 *
 *    Impact of losing samples?
 *    - sums in accumulators do not match total
 *      -> I'm thinking this is fine =S  don't really use the accumulator sums for anything except calculating their stats
 *    - accumulators are slightly inaccurate (1 out of every 2*wavelen samples lost, possible data to find wave)
 *      and slightly non-sensitive
 *
 * okay, given that i'm losing these tail samples, how do I properly include head and tail?
 * i calculate the shortest among the body.  do i need to includ head and tail?
 * -> head has an end, so can mess up samples
 *  -> tail will work fine; its end will track
 * but if head's end is too short, need to accomodate
 */
