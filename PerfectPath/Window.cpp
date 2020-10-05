#include <GL/glut.h>
#include <glm/glm.hpp>
#include <stdlib.h> 
#include <iostream> 
#include <algorithm> 

#include "Tracer.h"

const int width = 640;
const int height = 480;
glm::vec3** output;

void idle()
{
	glutPostRedisplay();
}

void display(void)
{

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, double(width), 0.0, double(height));

	glPointSize(1);

	glBegin(GL_POINTS); // begin draw
	if (output) {
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				glColor3f(std::min((int)output[col][row].x, 255), std::min((int)output[col][row].y, 255), std::min((int)output[col][row].z, 255));
				glVertex2i(col, row);
			}
		}
	}
	glEnd(); // end draw

	glFlush();
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE);
	glutInitWindowSize(width, height);
	glutCreateWindow("Perfect Path");
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
	output = render(640, 480);
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutMainLoop();
}
