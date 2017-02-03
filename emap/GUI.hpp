#pragma once

#include <memory>

#include <itpp/base/mat.h>
#include <itpp/base/vec.h>

class GUIWindow
{
public:
	// each vector will be charted as a line graph
	virtual void setLines(std::vector<itpp::vec> const &) = 0;

	// provided matrix will form a background image
	virtual void setImage(itpp::mat const &) = 0;

	// a single row is added to the bottom of the background image, scrolling the rest up
	virtual void addRow(itpp::vec const &) = 0;

	virtual std::pair<int, int> size() = 0;

	static std::unique_ptr<GUIWindow> create(unsigned width = 0, unsigned height = 0);

	virtual ~GUIWindow() = 0;
};
