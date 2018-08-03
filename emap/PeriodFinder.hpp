#pragma once

#include "GUI.hpp"
#include "Source.hpp"

/*
 * Precisely finding a period is important for extracting weak repeating signals from noise by recording for a long time.
 *
 * I'm finding there may be different stages to finding a precise average period:
 *
 * 1. Finding the period among the entire conceptual bandwidth.
 *        This moves from a concept of say "it happens around 4 times a second" to "it happens at 4.0 Hz"
 *
 *        Approaches include manually guessing a number, or looking for spectral peaks taken from a long or highly decimated recording.
 *
 * 2. Finding the period to the precision of the recording samplerate
 *        Moves from the concept of "4.0 Hz" to "The period is on average 600179 samples long".
 *
 * 3. Finding the period to subsample precision, to the maximum precision given by the number of periods recorded.
 *        Moves from "600179 samples long" to "There were exactly 4002 periods in this 1000 second recording."
 *
 *        One approach may involve solving this problem for a short recording time, guessing that the period is stable and summarizing the
 *        data without storing it all, then using the next block of data to refine what is already known.
 *
 * 4. Once the average period is known, the variance of the period may be explored.
 *        It is expected to follow somewhat predictable changes that could be graphed.
 *
 * I used to imagine something like recording just enough such that the signal became visible, then perhaps recording another block, and seeing
 * how the period alignment matched.  Once the expected form of the signal is known, these blocks can be slid around or combined in different
 * ways to approximate a curve showing how the period timing is changing, and where the start of each period is likely to be, even though a
 * single period is too quiet to be visible alone.
 */

class PeriodFinder : public Destination
{
public:
	PeriodFinder(double minFrequency, double maxFrequency, double precisionFrequency, Source & source, double tunedHertz);
	~PeriodFinder();
	void receiveQuadrature(cvec const & data, double samplingHertz, double tunedHertz, double dBGain, double unixSecondsCompleted, class Source & source);

private:
	void init(double rate);

	std::unique_ptr<GUIWindow> window;

	std::pair<double,double> rangeHz;
  double precisionHz;
	Source & source;
	double tunedHz;

	double expectedSamplingHz;
	double longestPeriod, shortestPeriod;

	cvec data;
  typename cvec::value_type dataSum;

	std::vector<cvec> values;
	std::vector<vec> weights;
	std::vector<int> offsets;
};
