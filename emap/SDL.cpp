#include "SDL.hpp"

#include "Main.hpp"

#include <algorithm>
#include <future>
#include <system_error>
#include <thread>

#define SDL_MAIN_HANDLED
#include <SDL.h>
#include <SDL_ttf.h>

static TTF_Font * font;

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
		TEXTURE,
		SCROLL,
		TEXT
	} subtype;
	SDLWindow * window;
	std::promise<void> * promise;

	void * dataPtr;

	SDL_CustomWindowEvent(SubType type, SDLWindow * window, void const * dataPtr = 0)
	: SDL_CommonEvent({sdlErr(TYPE)}), subtype(type), window(sdlErr(window)), promise(new std::promise<void>()), dataPtr(const_cast<void*>(dataPtr))
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


SDLWindow::SDLWindow(std::string title, int width, int height)
: textureData(0, SDL_DestroyTexture), textData(0, SDL_DestroyTexture), window(0, SDL_DestroyWindow), renderer(0, SDL_DestroyRenderer)
{
	auto dims = std::make_tuple(title, width, height);
	SDL_CustomWindowEvent(SDL_CustomWindowEvent::CREATE, this, &dims);
}

void SDLWindow::performCreate(std::tuple<std::string,int,int> const * dims)
{
	std::string title = std::get<0>(*dims);
	int width = std::get<1>(*dims);
	int height = std::get<2>(*dims);
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
		SDL_Window * window;
		SDL_Renderer * renderer;
		sdlErr( SDL_CreateWindowAndRenderer(width, height, flags, &window, &renderer) );
		SDL_SetWindowTitle(window, title.c_str());
		this->window.reset(window);
		this->renderer.reset(renderer);
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
		windows.erase(SDL_GetWindowID(window.get()));
	}
	textData.release();
	textureData.release();
	renderer.release();
	window.release();
}

std::pair<int, int> SDLWindow::size()
{
	std::pair<int, int> dims;
	SDL_CustomWindowEvent event(SDL_CustomWindowEvent::GET_SIZE, this, &dims);
	return dims;
}

void SDLWindow::performSize(std::pair<int, int> * dims)
{
	sdlErr( SDL_GetRendererOutputSize(renderer.get(), &dims->first, &dims->second) );
}

void SDLWindow::setLines(std::vector<itpp::vec> const & values)
{
	{
		std::lock_guard<std::mutex> dataLk(dataMtx);
		vectorData = values;
	}
	SDL_Event event;
	event.window.type = SDL_WINDOWEVENT;
	event.window.windowID = SDL_GetWindowID(window.get());
	event.window.event = SDL_WINDOWEVENT_EXPOSED;
	sdlErr( SDL_PushEvent(&event) );
}

void SDLWindow::setImage(itpp::mat const & values)
{
	SDL_CustomWindowEvent event(SDL_CustomWindowEvent::TEXTURE, this, &values);
}

void SDLWindow::addRow(itpp::vec const & values)
{
	SDL_CustomWindowEvent event(SDL_CustomWindowEvent::SCROLL, this, &values);
}

void SDLWindow::setText(std::string const & value)
{
	SDL_CustomWindowEvent event(SDL_CustomWindowEvent::TEXT, this, &value);
}

void SDLWindow::performScroll(itpp::vec const * values)
{
	std::pair<int, int> dims;
	if (!textureData) {
		performSize(&dims);
		textureData.reset( sdlErr( SDL_CreateTexture(renderer.get(), SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, dims.first, dims.second) )  );
	} else {
		sdlErr( SDL_QueryTexture(textureData.get(), 0, 0, &dims.first, &dims.second) );
	}

	uint8_t * pixels;
	int pitch;

	sdlErr( SDL_LockTexture(textureData.get(), 0, reinterpret_cast<void**>(&pixels), &pitch) );
	unsigned offset = pitch * (dims.second - 1);
	memmove(pixels, pixels + pitch, offset);
	uint8_t * pixel = pixels + offset;
	for (int x = 0; x < dims.first; ++ x) {
		int x2 = x * values->size() / dims.first;
		uint8_t c = values->_elem(x2) * (256.0 - std::numeric_limits<double>::epsilon() * 256.0);
		*(pixel ++) = c;
		*(pixel ++) = c;
		*(pixel ++) = c;
	}

	SDL_UnlockTexture(textureData.get());
	performDraw();
}

void SDLWindow::performTexture(itpp::mat const * values)
{
	int w = 0, h = 0;
	if (textureData)
		sdlErr( SDL_QueryTexture(textureData.get(), 0, 0, &w, &h) );
	int height = values->rows();
	if (w != values->cols() || h != height) {
		w = values->cols();
		h = height;
		textureData.release();
		textureData.reset(sdlErr( SDL_CreateTexture(renderer.get(), SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, w, h) ));
	}

	uint8_t * pixels;
	int pitch;

	sdlErr( SDL_LockTexture(textureData.get(), 0, reinterpret_cast<void**>(&pixels), &pitch) );

	for (int y = 0, i = 0; y < values->rows(); ++ y, pixels += pitch) {
		uint8_t * pixel = pixels;
		for (int x = 0; x < values->cols(); ++ x, ++ i) {
			uint8_t c = values->_elem(y, x) * (256.0 - std::numeric_limits<double>::epsilon() * 256.0);
			*(pixel ++) = c;
			*(pixel ++) = c;
			*(pixel ++) = c;
		}
	}

	SDL_UnlockTexture(textureData.get());
	performDraw();
}

void SDLWindow::performText(std::string const * string)
{
	static SDL_Color fg = SDL_Color{255,0,255,0};
	static SDL_Color bg = SDL_Color{0,0,0,0};

	std::unique_ptr<SDL_Surface,void (*)(SDL_Surface*)> surface(TTF_RenderUTF8_Shaded(font, string->c_str(), fg, bg), SDL_FreeSurface);

	textData.reset(sdlErr( SDL_CreateTextureFromSurface(renderer.get(), surface.get()) ));

	performDraw();
}

void SDLWindow::performDraw()
{
	std::lock_guard<std::mutex> dataLk(dataMtx);
	if (textureData) {
		sdlErr( SDL_RenderCopy(renderer.get(), textureData.get(), 0, 0) );
	} else {
		sdlErr( SDL_SetRenderDrawColor(renderer.get(), 0, 0, 0, SDL_ALPHA_OPAQUE) );
		sdlErr( SDL_RenderClear(renderer.get()) );
	}
	if (vectorData.size()) {
		for (auto & data : vectorData) {
			static std::vector<SDL_Point> points;
			int w, h;
			sdlErr( SDL_GetRendererOutputSize(renderer.get(), &w, &h) );
		
			points.resize(data.size());
			for (size_t i = 0; i < points.size(); ++ i) {
				points[i].x = w * i / (points.size() - 1);
				points[i].y = h * (data[i] + 1) / 2;
			}
		
			sdlErr( SDL_SetRenderDrawColor(renderer.get(), 255, 0, 255, SDL_ALPHA_OPAQUE) );
			sdlErr( SDL_RenderDrawLines(renderer.get(), &points[0], points.size()) );
		}
	}
	if (textData) {
		SDL_Rect rect;
		rect.x = rect.y = 0;
		sdlErr( SDL_QueryTexture(textData.get(), 0, 0, &rect.w, &rect.h) );
		sdlErr( SDL_RenderCopy(renderer.get(), textData.get(), &rect, &rect) );
	}
	SDL_RenderPresent(renderer.get());
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
						performCreate(reinterpret_cast<std::tuple<std::string,int,int> const *>(e.dataPtr));
						break;
					case SDL_CustomWindowEvent::DESTROY:
						performDestroy();
						break;
					case SDL_CustomWindowEvent::TEXTURE:
						performTexture(reinterpret_cast<itpp::mat const *>(e.dataPtr));
						break;
					case SDL_CustomWindowEvent::SCROLL:
						performScroll(reinterpret_cast<itpp::vec const *>(e.dataPtr));
						break;
					case SDL_CustomWindowEvent::TEXT:
						performText(reinterpret_cast<std::string const *>(e.dataPtr));
						break;
					case SDL_CustomWindowEvent::GET_SIZE:
						performSize(reinterpret_cast<std::pair<int,int> *>(e.dataPtr));
						break;
				}
				e.done();
			}
			break;
	}
}

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
		sdlErr( TTF_Init() );
		//sdlErr( SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1) );
		//sdlErr( SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 2) );

		SDL_CustomWindowEvent::TYPE = sdlErr( SDL_RegisterEvents(1) );

		font = sdlErr( TTF_OpenFont("LiberationMono-Bold.ttf", 12) );

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

		TTF_CloseFont(font);
		TTF_Quit();
		SDL_Quit();
	}

	std::thread th;
	std::promise<void> ready;
} sdlInstance;
