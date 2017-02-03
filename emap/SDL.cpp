#include "SDL.hpp"

// TODO: ensure all rendering happens in the SDL thread

#include "Main.hpp"

#include <algorithm>
#include <future>
#include <system_error>
#include <thread>

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

static uint32_t sdlErr(uint32_t r)
{
	if (r == 0xFFFFFFFF)
		sdlErr(-1);
	return r;
}


std::mutex windowsMtx;
std::unordered_map<uint32_t, SDLWindow *> windows;

struct SDL_CustomWindowEvent : public SDL_CommonEvent
{
	static uint32_t TYPE;
	enum SubType {
		CREATE,
		DESTROY,
		GET_SIZE,
		TEXTURE
	} subtype;
	SDLWindow * window;
	int w, h;
	std::promise<void> * promise;

	itpp::mat const * values;

	SDL_CustomWindowEvent(SubType type, SDLWindow * window, int w = 0, int h = 0, itpp::mat const * values = 0)
	: SDL_CommonEvent({sdlErr(TYPE)}), subtype(type), window(sdlErr(window)), w(w), h(h), promise(new std::promise<void>()), values(values)
	{
		auto future = promise->get_future();
		sdlErr( SDL_PushEvent(reinterpret_cast<SDL_Event*>(this)) );
		future.wait();
		delete promise;
		promise = 0;
	}

	SDL_Event * event()
	{
		return reinterpret_cast<SDL_Event*>(this);
	}

	void done()
	{
		std::promise<void> promise = std::move(*this->promise);
		promise.set_value();
	}

	void except()
	{
		std::promise<void> promise = std::move(*this->promise);
		promise.set_exception(std::current_exception());
	}
};

uint32_t SDL_CustomWindowEvent::TYPE = -1;


SDLWindow::SDLWindow(unsigned width, unsigned height)
: textureData(0)
{
	SDL_CustomWindowEvent(SDL_CustomWindowEvent::CREATE, this, width, height);
}

void SDLWindow::performCreate(int width, int height)
{
	int flags = 0;
	if (width == 0 || height == 0) {
		SDL_Rect bounds;
		sdlErr( SDL_GetDisplayBounds(0, &bounds) );
		width = bounds.w / 2;
		height = bounds.h / 2;
		flags |= SDL_WINDOW_RESIZABLE;
	}
	{
		std::lock_guard<std::mutex> windowsLk(windowsMtx);
		sdlErr( SDL_CreateWindowAndRenderer(width, height, flags, &window, &renderer) );
		windows[SDL_GetWindowID(window)] = this;
	}
}

SDLWindow::~SDLWindow()
{
	SDL_CustomWindowEvent(SDL_CustomWindowEvent::DESTROY, this);
}

void SDLWindow::performDestroy()
{
	{
		std::lock_guard<std::mutex> windowsLk(windowsMtx);
		windows.erase(SDL_GetWindowID(window));
	}
	if (textureData) {
		SDL_DestroyTexture(textureData);
		textureData = 0;
	}
	SDL_DestroyRenderer(renderer);
	renderer = 0;
	SDL_DestroyWindow(window);
	window = 0;
}

std::pair<unsigned, unsigned> SDLWindow::size()
{
	SDL_CustomWindowEvent event(SDL_CustomWindowEvent::GET_SIZE, this);
	return std::make_pair(event.w, event.h);
}

void SDLWindow::performSize(int & w, int & h)
{
	sdlErr( SDL_GetRendererOutputSize(renderer, &w, &h) );
}

void SDLWindow::draw(std::vector<itpp::vec> const & values)
{
	{
		std::lock_guard<std::mutex> dataLk(dataMtx);
		vectorData = values;
	}
	SDL_Event event;
	event.window.type = SDL_WINDOWEVENT;
	event.window.windowID = SDL_GetWindowID(window);
	event.window.event = SDL_WINDOWEVENT_EXPOSED;
	sdlErr( SDL_PushEvent(&event) );
}

void SDLWindow::draw(itpp::mat const & values)
{
	SDL_CustomWindowEvent event(SDL_CustomWindowEvent::TEXTURE, this, 0, 0, &values);
}

void SDLWindow::performTexture(itpp::mat const * values)
{
	int w = 0, h = 0;
	if (textureData)
		sdlErr( SDL_QueryTexture(textureData, 0, 0, &w, &h) );
	if (w != values->rows() || h != values->cols()) {
		w = values->cols();
		h = values->rows();
		SDL_DestroyTexture(textureData);
		textureData = 0;
		textureData = sdlErr( SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, w, h) );
	}

	uint8_t * pixels;
	int pitch;

	sdlErr( SDL_LockTexture(textureData, 0, reinterpret_cast<void**>(&pixels), &pitch) );

	for (int y = 0; y < values->rows(); ++ y, pixels += pitch) {
		uint8_t * pixel = pixels;
		for (int x = 0; x < values->cols(); ++ x) {
			uint8_t c = std::max(std::min(values->get(y, x) * 256.0, 255.0), 0.0);
			*(pixel ++) = c;
			*(pixel ++) = c;
			*(pixel ++) = c;
		}
	}

	SDL_UnlockTexture(textureData);
	performDraw();
}

void SDLWindow::performDraw()
{
	std::lock_guard<std::mutex> dataLk(dataMtx);
	if (textureData) {
		sdlErr( SDL_RenderCopy(renderer, textureData, 0, 0) );
	} else {
		sdlErr( SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE) );
		sdlErr( SDL_RenderClear(renderer) );
	}
	if (vectorData.size()) {
		for (auto & data : vectorData) {
			static std::vector<SDL_Point> points;
			int w, h;
			sdlErr( SDL_GetRendererOutputSize(renderer, &w, &h) );
		
			points.resize(data.size());
			for (size_t i = 0; i < points.size(); ++ i) {
				points[i].x = w * i / (points.size() - 1);
				points[i].y = h * (data[i] + 1) / 2;
			}
		
			sdlErr( SDL_SetRenderDrawColor(renderer, 255, 0, 255, SDL_ALPHA_OPAQUE) );
			sdlErr( SDL_RenderDrawLines(renderer, &points[0], points.size()) );
		}
	}
	SDL_RenderPresent(renderer);
}

void SDLWindow::receiveSDLEvent(SDL_Event & event)
{
	switch (event.type)
	{
		case SDL_WINDOWEVENT:
			switch (event.window.event)
			{
			//case SDL_WINDOWEVENT_SIZE_CHANGED:
			case SDL_WINDOWEVENT_EXPOSED:
				performDraw();
				break;
			case SDL_WINDOWEVENT_CLOSE:
				{
					SDL_Event e;
					e.type = SDL_QUIT;
					sdlErr( SDL_PushEvent(&e) );
				}
				break;
			}
			break;
		default: if (event.type == SDL_CustomWindowEvent::TYPE)
			{
				SDL_CustomWindowEvent & e = *reinterpret_cast<SDL_CustomWindowEvent *>(& event);
				switch(e.subtype) {
					case SDL_CustomWindowEvent::CREATE:
						performCreate(e.w, e.h);
						break;
					case SDL_CustomWindowEvent::DESTROY:
						performDestroy();
						break;
					case SDL_CustomWindowEvent::TEXTURE:
						performTexture(e.values);
						break;
					case SDL_CustomWindowEvent::GET_SIZE:
						performSize(e.w, e.h);
						break;
				}
				e.done();
			}
			break;
	}
}

#include <iostream>

static class SDL
{
public:
	SDL()
	{
		th = std::thread(&SDL::run, this);
		ready.get_future().wait();
	}

	~SDL()
	{
		SDL_Event event;
		event.type = SDL_QUIT;
		sdlErr( SDL_PushEvent(&event) );
		th.join();
	}
private:
	void run()
	{
		SDL_SetMainReady();
		sdlErr( SDL_Init(SDL_INIT_VIDEO) );
		//sdlErr( SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1) );
		//sdlErr( SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 2) );

		SDL_CustomWindowEvent::TYPE = sdlErr( SDL_RegisterEvents(1) );

		ready.set_value();

		SDL_Event event;
		bool quit = false;
		bool quitting = false;
		do {
			sdlErr( 0 != SDL_WaitEvent(&event) );
			switch (event.type)
			{
			case SDL_WINDOWEVENT:
				{
					if (quitting)
						break;
					SDLWindow * window;
					{
						std::lock_guard<std::mutex> windowsLk(windowsMtx);
						auto windowIt = windows.find(event.window.windowID);
						if (windowIt == windows.end())
							break;
						window = windowIt->second;
					}
					window->receiveSDLEvent(event);
				}
				break;
			case SDL_QUIT:
				{
					std::lock_guard<std::mutex> windowsLk(windowsMtx);
					quit = windows.empty();
					if (quitting)
						break;
					Main::instance().stop();
					quitting = true;
				}
				break;
			default: if (event.type == SDL_CustomWindowEvent::TYPE)
				{
					SDL_CustomWindowEvent * e = reinterpret_cast<SDL_CustomWindowEvent *>(&event);
					try {
						if (quitting && e->subtype != SDL_CustomWindowEvent::DESTROY)
							throw std::invalid_argument("quitting");
						SDLWindow * window = e->window;
						window->receiveSDLEvent(event);
					} catch (...) {
						e->except();
					}
					break;
				}
			}
		} while (!quit);

		SDL_Quit();
	}

	std::thread th;
	std::promise<void> ready;
} sdlInstance;
