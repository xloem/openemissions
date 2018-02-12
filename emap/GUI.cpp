#include "GUI.hpp"

#include "SDL.hpp"

std::unique_ptr<GUIWindow> GUIWindow::create(std::string title, unsigned width, unsigned height)
{
	return std::unique_ptr<GUIWindow>(new SDLWindow(title, width, height));
}

GUIWindow::~GUIWindow()
{ }
