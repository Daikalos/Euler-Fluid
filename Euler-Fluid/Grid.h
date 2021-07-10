#pragma once

#include <SFML/Graphics.hpp>
#include <unordered_map>
#include "VecUtil.h"
#include "Cell.h"

template<typename T> class Grid
{
public:
	Grid(int grid_left, int grid_top, int grid_right, int grid_bot, int cont_width, int cont_height);
	~Grid();

	inline Container<T>* at_pos(const sf::Vector2f& position) const
	{
		const sf::Vector2i& pos = ((sf::Vector2i)position - gridRect.top_left) / contDims;

		if (!within_grid(pos))
			return nullptr;

		return at_pos(pos);
	}
	inline Container<T>* at_pos(const T& item) const
	{
		const Boid* boid = dynamic_cast<const Boid*>(&item);

		if (boid == nullptr)
			return nullptr;

		const sf::Vector2i& pos = ((sf::Vector2i)boid->get_origin() - gridRect.top_left) / contDims;

		if (!within_grid(pos))
			return nullptr;

		return at_pos(pos);
	}

private:
	inline Container<T>* at_pos(int x, int y) const
	{
		return &containers[x + y * width];
	}
	inline Container<T>* at_pos(const sf::Vector2i& position) const
	{
		return &containers[position.x + position.y * width];
	}

	inline bool within_grid(const sf::Vector2i& pos) const
	{
		return !(pos.x < 0 || pos.y < 0 || pos.x >= width || pos.y >= height);
	}

private:
	Container<T>* containers;
	sf::Vector2i contDims;

	std::unordered_map<const T*, Container<T>*> items;

	unsigned short width, height;
	Rect_i gridRect;

private:
	Grid() = delete;
	Grid(const Grid& rhs) = delete;
	Grid& operator=(const Grid& rhs) = delete;
};

