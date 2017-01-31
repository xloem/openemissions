#pragma once

#include <memory>

#include <itpp/base/mat.h>
#include <itpp/base/vec.h>

class GUIWindow
{
public:
	virtual void draw(itpp::vec const &) = 0;
	virtual void draw(itpp::mat const &) = 0;

	static std::unique_ptr<GUIWindow> create(unsigned width = 0, unsigned height = 0);

	virtual ~GUIWindow() = 0;
};
