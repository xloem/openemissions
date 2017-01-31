#include "GUI.hpp"

#include "SDL.hpp"

std::unique_ptr<GUIWindow> GUIWindow::create(unsigned width, unsigned height)
{
	return std::unique_ptr<GUIWindow>(new SDLWindow(width, height));
}

GUIWindow::~GUIWindow()
{ }
