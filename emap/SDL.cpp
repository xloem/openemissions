#include "SDL.hpp"

#include <algorithm>
#include <system_error>

#define SDL_MAIN_HANDLED
#include <SDL.h>

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
: textureData(0)
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
	renderer = 0;
	SDL_DestroyWindow(window);
	window = 0;
	if (textureData) {
		SDL_DestroyTexture(textureData);
		textureData = 0;
	}
}

std::pair<unsigned, unsigned> SDLWindow::size()
{
	int w, h;
	sdlErr( SDL_GetRendererOutputSize(renderer, &w, &h) );
	return std::make_pair<unsigned, unsigned>(w, h);
}

void SDLWindow::draw(itpp::vec const & values)
{
	vectorData = values;
	draw();
}

void SDLWindow::draw(itpp::mat const & values)
{
	int w = 0, h = 0;
	if (textureData)
		sdlErr( SDL_QueryTexture(textureData, 0, 0, &w, &h) );
	if (w != values.cols() || h != values.rows()) {
		w = values.cols();
		h = values.rows();
		SDL_DestroyTexture(textureData);
		textureData = 0;
		textureData = sdlErr( SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, w, h) );
	}

	void * pixels;
	int pitch;
	sdlErr( SDL_LockTexture(textureData, 0, &pixels, &pitch) );
	for (int y = 0; y < h; ++ y, pixels = static_cast<uint8_t *>(pixels) + pitch) {
		uint8_t * pixel = static_cast<uint8_t *>(pixels);
		for (int x = 0; x < w; ++ x) {
			uint8_t c = std::max(std::min(values(y, x) * 256.0, 255.0), 0.0);
			*(pixel ++) = c;
			*(pixel ++) = c;
			*(pixel ++) = c;
		}
	}
	SDL_UnlockTexture(textureData);
	draw();
}

void SDLWindow::draw() {
	if (textureData) {
		sdlErr( SDL_RenderCopy(renderer, textureData, 0, 0) );
	} else {
		sdlErr( SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE) );
		sdlErr( SDL_RenderClear(renderer) );
	}
	if (vectorData.size()) {
		static std::vector<SDL_Point> points;
		auto dims = this->size();
	
		points.resize(vectorData.size());
		for (size_t i = 0; i < points.size(); ++ i) {
			points[i].x = dims.first * i / (points.size() - 1);
			points[i].y = dims.second * vectorData[i];
		}
	
		sdlErr( SDL_SetRenderDrawColor(renderer, 255, 0, 255, SDL_ALPHA_OPAQUE) );
		sdlErr( SDL_RenderDrawLines(renderer, &points[0], points.size()) );
	}
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
