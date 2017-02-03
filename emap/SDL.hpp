#pragma once

#include "GUI.hpp"

#include <mutex>
#include <unordered_map>
#include <itpp/base/mat.h>
#include <itpp/base/vec.h>

class SDLWindow : public GUIWindow
{
friend class SDL;
public:
	SDLWindow(unsigned width = 0, unsigned height = 0);
	~SDLWindow();

	void draw(std::vector<itpp::vec> const &);
	void draw(itpp::mat const &);
	std::pair<unsigned, unsigned> size();

private:
	std::mutex dataMtx;
	std::vector<itpp::vec> vectorData;
	struct SDL_Texture * textureData;

	void performCreate(int width, int height);
	void performSize(int & width, int & height);
	void performTexture(itpp::mat const * values);
	void performDraw();
	void performDestroy();

	struct SDL_Window * window;
	struct SDL_Renderer * renderer;

	void receiveSDLEvent(union SDL_Event & event);
};
