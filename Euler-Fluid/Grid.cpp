#include "Grid.h"

template<typename T>
Grid<T>::Grid(int grid_left, int grid_top, int grid_right, int grid_bot, int cont_width, int cont_height)
	: contDims({ cont_width, cont_height }), gridRect(grid_left, grid_top, grid_right, grid_bot)
{
	width = gridRect.width() / cont_width;
    height = gridRect.height() / cont_height;

	containers = new Container<T>[(long)width * height];

	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			containers[x + y * (long)width] = Container<T>(Rect_i(
				 x * cont_width + grid_left,			  y * cont_height + grid_top,
				 x * cont_width + cont_width + grid_left, y * cont_height + cont_height + grid_top));
		}
	}
}

template<typename T>
Grid<T>::~Grid()
{
	delete[] containers;
}
