#pragma once

#include <SFML/Graphics.hpp>

template<typename T> struct Rect
{
	Rect()
	{
		top_left = bot_right = top_right = bot_left = {};
	}
	Rect(T left, T top, T right, T bot)
		: top_left({ left, top }), bot_right({ right, bot })
	{
		top_right = { bot_right.x, top_left.y };
		bot_left = { top_left.x, bot_right.y };
	};
	Rect(sf::Vector2<T> top_left, sf::Vector2<T> bot_right)
		: top_left(top_left), bot_right(bot_right)
	{
		top_right = { bot_right.x, top_left.y };
		bot_left = { top_left.x, bot_right.y };
	};

	inline T width() const { return (right - left); }
	inline T height() const { return (bot - top); }

	union
	{
		sf::Vector2<T> top_left;
		struct { T left, top; };
	};
	union
	{
		sf::Vector2<T> bot_right;
		struct { T right, bot; };
	};

	sf::Vector2<T>
		top_right,
		bot_left;
};

typedef Rect<int> Rect_i;