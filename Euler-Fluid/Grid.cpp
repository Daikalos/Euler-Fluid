#include "Grid.h"

Grid::Grid(int grid_left, int grid_top, int grid_right, int grid_bot, int cont_width, int cont_height)
	: cellDims({ cont_width, cont_height }), gridRect(grid_left, grid_top, grid_right, grid_bot)
{
	width = gridRect.width() / cont_width;
    height = gridRect.height() / cont_height;

	count = width * height;

	cells = new Cell[count];

	for (int y = 0; y < height; ++y)
		for (int x = 0; x < width; ++x)
		{
			cells[x + y * (long)width] = Cell();
		}
}

Grid::~Grid()
{
	delete[] cells;
}

void Grid::update(float deltaTime)
{
	std::for_each(
		std::execution::par_unseq,
		cells,
		cells + count,
		[&](Cell& cell)
		{
			cell.update(deltaTime);
		}
	);
}
