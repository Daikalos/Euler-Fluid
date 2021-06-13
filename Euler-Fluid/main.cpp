#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <gl/GLU.h>

#include "Camera.h"
#include "InputHandler.h"

int main()
{
	sf::RenderWindow window(sf::VideoMode(2240, 1260), "Euler Fluid");
	window.setActive(true);

	sf::CircleShape circle(100.0f);
	circle.setFillColor(sf::Color::Red);
	Camera camera(window);
	InputHandler inputHandler;

	sf::Clock clock;
	float deltaTime = FLT_EPSILON;

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glScalef(1.0f, -1.0f, 1.0f);
	gluOrtho2D(0, window.getSize().x, 0, window.getSize().y);

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
				gluOrtho2D(0, window.getSize().x, 0, window.getSize().y);
				glMatrixMode(GL_MODELVIEW);
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