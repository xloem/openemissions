#include "FLAC_XMP.hpp"

#include <cstdlib>

#include "FLAC++/decoder.h"
#include "tinyxml2.h"

class FLACFileDecoder : public FLAC::Decoder::File
{
public:
	FLACFileDecoder(std::string file, FLAC_XMP & source)
	: source(source)
	{
		init(file);
	}
	~FLACFileDecoder()
	{
		finish();
	}

protected:
	FLAC__StreamDecoderWriteStatus write_callback(const FLAC__Frame *frame, const FLAC__int32 * const buffer[])
	{
		if (frame->header.channels != 2)
			return FLAC__STREAM_DECODER_WRITE_STATUS_ABORT;

		unsigned samples = frame->header.blocksize;
		itpp::cvec data(samples);

		for (unsigned i = 0; i < samples; ++ i) {
			data[i].real(buffer[0][i] / 127.0);
			data[i].imag(buffer[1][i] / 127.0);
		}

		source.processSamples(data);
		
		return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
	}

	void error_callback(FLAC__StreamDecoderErrorStatus status)
	{
	}

private:
	FLAC_XMP & source;
};

#include <iostream>

FLAC_XMP::FLAC_XMP(std::string flacPath, std::string xmpPath)
: flacDecoder(new FLACFileDecoder(flacPath, *this)),
  xmlDocument(new tinyxml2::XMLDocument()),
  sampleRate(0), tuneFrequency(0), resultTime(0), xmpStart(0), xmpDuration(0), flacOffset(0)
{
	xmlDocument->LoadFile(xmpPath.c_str());
	auto xmpmeta = xmlDocument->FirstChildElement("x:xmpmeta");
	auto RDF = xmpmeta->FirstChildElement("rdf:RDF");
	auto descr = RDF->FirstChildElement("rdf:Description");
	auto tracks = descr->FirstChildElement("xmpDM:Tracks");
	auto bag = tracks->FirstChildElement("rdf:Bag");
	auto resource = bag->FirstChildElement("rdf:li");
	nextElement = resource->FirstChildElement();

	flacDecoder->process_until_end_of_stream();
}

FLAC_XMP::~FLAC_XMP()
{
}

void FLAC_XMP::processSamples(itpp::cvec const & data)
{
	size_t flacIndex = 0;

	while (flacOffset >= xmpStart + xmpDuration) {
		// process xmp
		if (std::string("xmpDM:frameRate") == nextElement->Name()) {
			char const * frameRateText = nextElement->GetText();
			while (*frameRateText && (*frameRateText < '0' || *frameRateText > '9'))
				++ frameRateText;
			sampleRate = std::strtod(frameRateText, 0);
		} else if (std::string("xmpDM:startTime") == nextElement->Name()) {
			auto root = nextElement->Parent()->ToElement();
			root->FirstChildElement("xmpDM:startTime")->QueryUnsignedText(&xmpStart);
			root->FirstChildElement("xmpDM:duration")->QueryUnsignedText(&xmpDuration);
			auto seq = root->FirstChildElement("xmpDM:cuePointParams")->FirstChildElement("rdf:Seq");
			for (auto elem = seq->FirstChildElement(); elem; elem = elem->NextSiblingElement()) {
				if (std::string("ResultTime") == elem->FirstChildElement("xmpDM:key")->GetText()) {
					elem->FirstChildElement("xmpDM:value")->QueryDoubleText(&resultTime);
				} else if (std::string("TuneFrequency") == elem->FirstChildElement("xmpDM:key")->GetText()) {
					elem->FirstChildElement("xmpDM:value")->QueryDoubleText(&tuneFrequency);
				}
			}
			std::cout << "xmpStart: " << xmpStart << " xmpDuration: " << xmpDuration << " sampleRate: " << sampleRate << " frequency: " << tuneFrequency << std::endl;
		}
		tinyxml2::XMLElement * next = nextElement->FirstChildElement();
		if (!next)
			next = nextElement->NextSiblingElement();
		if (!next)
			next = nextElement->Parent()->ToElement();
		nextElement = next;
	}
	
	unsigned count = xmpDuration;
	if (flacIndex + count < static_cast<unsigned>(data.size()))
		count = data.size() - flacIndex;
	
	
}
