#pragma once

#include "Source.hpp"

#include "itpp/base/vec.h"

namespace tinyxml2 { class XMLDocument; class XMLElement; }

class FLAC_XMP : public Source
{
friend class FLACFileDecoder;
public:
	FLAC_XMP(std::string flacPath, std::string xmpPath);
	~FLAC_XMP();

private:
	void processSamples(itpp::cvec const & data);

	std::unique_ptr<class FLACFileDecoder> flacDecoder;
	std::unique_ptr<tinyxml2::XMLDocument> xmlDocument;

	tinyxml2::XMLElement * nextElement;

	double sampleRate;
	double tuneFrequency;
	double resultTime;
	unsigned xmpStart;
	unsigned xmpDuration;
	size_t flacOffset;
};
