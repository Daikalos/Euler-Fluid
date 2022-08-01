#pragma once

#include <SFML/Graphics.hpp>

#include <iostream>
#include <unordered_map>

#include "Config.h"

class ResourceManager
{
public:
	ResourceManager();
	~ResourceManager();

	sf::Texture* request_texture(std::string name, sf::Texture* fallback = nullptr) const;
	sf::Font* request_font(std::string name, sf::Font* fallback = nullptr) const;

	void load_texture(std::string name, std::string path);
	void load_font(std::string name, std::string path);

	void load_textures();
	void load_fonts();

	void clean_up();

private:
	std::unordered_map<std::string, std::shared_ptr<sf::Texture>> textures;
	std::unordered_map<std::string, std::shared_ptr<sf::Font>> fonts;
};

