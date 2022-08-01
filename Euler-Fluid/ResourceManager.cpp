#include "ResourceManager.h"

ResourceManager::ResourceManager()
{
}

ResourceManager::~ResourceManager()
{
	clean_up();
}

sf::Texture* ResourceManager::request_texture(std::string name, sf::Texture* fallback) const
{
	auto it = textures.find(name);

	if (it != textures.end())
		return it->second.get();
	
	return fallback;
}

sf::Font* ResourceManager::request_font(std::string name, sf::Font* fallback) const
{
	auto it = fonts.find(name);

	if (it != fonts.end())
		return it->second.get();

	return fallback;
}
void ResourceManager::load_textures()
{

}
void ResourceManager::load_fonts()
{

}

void ResourceManager::clean_up()
{
	textures.clear();
	fonts.clear();
}

void ResourceManager::load_texture(std::string name, std::string path)
{
	std::shared_ptr<sf::Texture> texture = std::make_shared<sf::Texture>();

	if (!texture->loadFromFile(path))
	{
		textures.erase(name);
		return;
	}

	textures[name] = texture;
}
void ResourceManager::load_font(std::string name, std::string path)
{
	std::shared_ptr<sf::Font> font = std::make_shared<sf::Font>();

	if (!font->loadFromFile(path))
	{
		fonts.erase(name);
		return;
	}

	fonts[name] = font;
}
