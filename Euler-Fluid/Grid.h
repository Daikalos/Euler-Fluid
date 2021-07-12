#pragma once

#include <SFML/Graphics.hpp>
#include <execution>

#include "VecUtil.h"
#include "Rectangle.h"
#include "Cell.h"

class Grid
{
public:
	Grid(int grid_left, int grid_top, int grid_right, int grid_bot, int cont_width, int cont_height);
	~Grid();

	inline Cell* at_pos(const sf::Vector2f& pos) const
	{
		const sf::Vector2i& position = ((sf::Vector2i)pos - gridRect.top_left) / cellDims;

		if (!within_grid(position))
			return nullptr;

		return at_pos(position);
	}

	void update(float deltaTime);

private:
	inline Cell* at_pos(int x, int y) const
	{
		return &cells[x + y * width];
	}
	inline Cell* at_pos(const sf::Vector2i& pos) const
	{
		return &cells[pos.x + pos.y * width];
	}

	inline bool within_grid(const sf::Vector2i& pos) const
	{
		return !(pos.x < 0 || pos.y < 0 || pos.x >= width || pos.y >= height);
	}

private:
	Cell* cells;
	sf::Vector2i cellDims;

	unsigned short width, height, count;
	Rect_i gridRect;

private:
	Grid() = delete;
	Grid(const Grid& rhs) = delete;
	Grid& operator=(const Grid& rhs) = delete;
};

