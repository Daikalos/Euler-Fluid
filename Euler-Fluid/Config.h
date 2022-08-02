#pragma once

#include <SFML/Graphics.hpp>
#include <nlohmann/json.hpp>
#include <fstream>

static const std::string FILE_NAME = "config";

enum Reconstruct
{
	RFluid
};

struct Config
{
	int scale;

public:
	Config();
	~Config();

	void load();
	std::vector<Reconstruct> refresh(Config& prev);

private:
	void load_var(nlohmann::json& json);

	sf::Vector3f convert(std::string strColor) const
	{
		std::stringstream stream(strColor);
		std::string segment;
		std::vector<std::string> values;

		while (std::getline(stream, segment, ' '))
			values.push_back(segment);

		sf::Vector3f color(1.0f, 1.0f, 1.0f);

		try
		{
			color.x = std::stof(values[0]) / 255.0f;
			color.y = std::stof(values[1]) / 255.0f;
			color.z = std::stof(values[2]) / 255.0f;
		}
		catch (std::exception e) {}

		return color;
	}
};
