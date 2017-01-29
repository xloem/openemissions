#include "SDL.hpp"

#define SDL_MAIN_HANDLED
#include <SDL.h>
#include <system_error>

static int sdlErr(int r)
{
	if (r < 0)
		throw std::system_error(r, std::generic_category(), SDL_GetError());
	return r;
}

static void sdlErr(bool r)
{
	if (!r)
		sdlErr(-1);
}

template <typename T>
static T * sdlErr(T * r)
{
	if (r == 0)
		sdlErr(-1);
	return r;
}

std::unordered_map<uint32_t, SDLWindow *> SDLWindow::windows;

SDLWindow::SDLWindow(int width, int height)
{
	int flags = 0;
	if (width == 0 || height == 0) {
		SDL_Rect bounds;
		sdlErr( SDL_GetDisplayBounds(0, &bounds) );
		width = bounds.w / 2;
		height = bounds.h / 2;
		flags |= SDL_WINDOW_RESIZABLE;
	}
	sdlErr( SDL_CreateWindowAndRenderer(width, height, flags, &window, &renderer) );
	windows[SDL_GetWindowID(window)] = this;
}
SDLWindow::~SDLWindow()
{
	windows.erase(SDL_GetWindowID(window));
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
}

std::pair<unsigned, unsigned> SDLWindow::size()
{
	int w, h;
	sdlErr( SDL_GetRendererOutputSize(renderer, &w, &h) );
	return std::make_pair<unsigned, unsigned>(w, h);
}

void SDLWindow::draw(std::vector<double> const & values)
{
	vectorData = values;
	draw();
}

void SDLWindow::draw() {
	static std::vector<SDL_Point> points;
	auto dims = this->size();

	points.resize(vectorData.size());
	for (size_t i = 0; i < points.size(); ++ i) {
		points[i].x = dims.first * i / (points.size() - 1);
		points[i].y = dims.second * vectorData[i];
	}

	sdlErr( SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE) );
	sdlErr( SDL_RenderClear(renderer) );
	sdlErr( SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE) );
	sdlErr( SDL_RenderDrawLines(renderer, &points[0], points.size()) );
	SDL_RenderPresent(renderer);
}

void SDLWindow::receiveSDLEvent(SDL_WindowEvent & event)
{
	switch (event.event)
	{
	case SDL_WINDOWEVENT_SIZE_CHANGED:
		draw();
		break;
	}
}

void SDLEvents::loop()
{
	SDL_Event event;
	for (;;) {
		sdlErr( 0 != SDL_WaitEvent(&event) );
		switch (event.type)
		{
		case SDL_QUIT:
			return;
		case SDL_WINDOWEVENT:
			SDLWindow::windows[event.window.windowID]->receiveSDLEvent(event.window);
			break;
		}
	}
}

static class SDL
{
public:
	SDL()
	{
		SDL_SetMainReady();
		sdlErr( SDL_Init(SDL_INIT_VIDEO) );
		//sdlErr( SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1) );
		//sdlErr( SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 2) );
	}
	~SDL()
	{
		SDL_Quit();
	}
} sdlInstance;
