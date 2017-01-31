#pragma once

#include <mutex>
#include <unordered_map>
#include <itpp/base/mat.h>
#include <itpp/base/vec.h>

class SDLWindow
{
friend class SDL;
public:
	SDLWindow(unsigned width = 0, unsigned height = 0);
	~SDLWindow();

	void draw(itpp::vec const &);
	void draw(itpp::mat const &);
	std::pair<unsigned, unsigned> size();

private:
	std::mutex dataMtx;
	itpp::vec vectorData;
	struct SDL_Texture * textureData;

	void performCreate(int width, int height);
	void performSize(int & width, int & height);
	void performLock(int width, int height, void * * pixels, int * pitch);
	void performUnlock();
	void performDraw();
	void performDestroy();

	struct SDL_Window * window;
	struct SDL_Renderer * renderer;

	void receiveSDLEvent(union SDL_Event & event);
};
