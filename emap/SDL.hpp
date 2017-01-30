#pragma once

#include <unordered_map>
#include <itpp/base/mat.h>
#include <itpp/base/vec.h>

class SDLWindow
{
friend class SDLEvents;
public:
	SDLWindow(int width = 0, int height = 0);
	~SDLWindow();

	void draw(itpp::vec const &);
	void draw(itpp::mat const &);
	std::pair<unsigned, unsigned> size();

private:
	itpp::vec vectorData;
	struct SDL_Texture * textureData;
	void draw();

	struct SDL_Window * window;
	struct SDL_Renderer * renderer;

	static std::unordered_map<uint32_t, SDLWindow *> windows;
	void receiveSDLEvent(struct SDL_WindowEvent & event);
};

class SDLEvents
{
public:
	static void loop();
};
