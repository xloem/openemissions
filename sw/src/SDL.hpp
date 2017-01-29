#pragma once

#include <vector>
#include <unordered_map>

class SDLWindow
{
friend class SDLEvents;
public:
	SDLWindow(int width = 0, int height = 0);
	~SDLWindow();

	void draw(std::vector<double> const &);
	std::pair<unsigned, unsigned> size();

private:
	std::vector<double> vectorData;
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
