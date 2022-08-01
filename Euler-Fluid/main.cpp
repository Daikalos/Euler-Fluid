#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>

#include "Camera.h"
#include "InputHandler.h"
#include "Grid.h"
#include "Cell.h"

int main()
{
	srand(time(NULL));

	std::vector<sf::VideoMode> fullscreen_modes = sf::VideoMode::getFullscreenModes();

	if (fullscreen_modes.size() == 0)
		return -1;

	sf::VideoMode video_mode = fullscreen_modes.front();

	if (!video_mode.isValid())
		return -1;

	sf::RenderWindow window(video_mode, "Euler Fluid"); //sf::Style::Fullscreen);

	if (!window.setActive(true))
		return -1;

	window.setVerticalSyncEnabled(true);

	sf::CircleShape circle(100.0f);
	circle.setFillColor(sf::Color::Red);

	Camera camera(window);
	InputHandler inputHandler;

	sf::Clock clock;
	float deltaTime = FLT_EPSILON;

	Grid grid(0, 0, window.getSize().x, window.getSize().y, 64, 64);

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glScalef(1.0f, -1.0f, 1.0f);
	glOrtho(0, window.getSize().x, 0, window.getSize().y, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	int i = 0;
	while (window.isOpen())
	{
		deltaTime = clock.restart().asSeconds();

		inputHandler.update();

		sf::Event event;
		while (window.pollEvent(event))
		{
			switch (event.type)
			{
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::Resized:
				glViewport(0, 0, window.getSize().x, window.getSize().y);
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				glScalef(1.0f, -1.0f, 1.0f);
				glOrtho(0, window.getSize().x, 0, window.getSize().y, -1.0, 1.0);
				glMatrixMode(GL_MODELVIEW);

				camera.set_position((sf::Vector2f)window.getSize() / 2.0f);

				break;
			case sf::Event::MouseWheelScrolled:
				inputHandler.set_scrollDelta(event.mouseWheelScroll.delta);
				break;
			}
		}


		camera.update(inputHandler);
		
		window.setView(camera.get_view());

		window.clear();

		window.draw(circle);

		window.display();
	}

	return 0;
}