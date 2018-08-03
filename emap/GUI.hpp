#pragma once

#include <memory>
#include <string>

#include "Common.hpp"

class GUIWindow
{
public:
	// each vector will be charted as a line graph
	virtual void setLines(std::vector<vec> const &) = 0;

	// provided matrix will form a background image
	virtual void setImage(mat const &) = 0;

	// a single row is added to the bottom of the background image, scrolling the rest up
	virtual void addRow(vec const &) = 0;

	// text may be displayed on top
	virtual void setText(std::string const &) = 0;

	virtual std::pair<int, int> size() = 0;

	static std::unique_ptr<GUIWindow> create(std::string title = "", unsigned width = 0, unsigned height = 0);

	virtual ~GUIWindow() = 0;
};
