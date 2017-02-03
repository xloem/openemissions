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
	SDLWindow(int width = 0, int height = 0);
	~SDLWindow();

	void setLines(std::vector<itpp::vec> const &);
	void setImage(itpp::mat const &);
	void addRow(itpp::vec const &);
	void setText(std::string const &);

	std::pair<int, int> size();

private:
	std::mutex dataMtx;
	std::vector<itpp::vec> vectorData;
	std::unique_ptr<struct SDL_Texture,void(*)(struct SDL_Texture*)> textureData;
	std::unique_ptr<struct SDL_Texture,void(*)(struct SDL_Texture*)> textData;

	void performCreate(std::pair<int,int> const * dims);
	void performSize(std::pair<int,int> * dims);
	void performTexture(itpp::mat const * values);
	void performScroll(itpp::vec const * values);
	void performText(std::string const * string);
	void performDraw();
	void performDestroy();

	std::unique_ptr<struct SDL_Window,void(*)(struct SDL_Window*)> window;
	std::unique_ptr<struct SDL_Renderer,void(*)(struct SDL_Renderer*)> renderer;

	void receiveSDLEvent(union SDL_Event & event);
};
