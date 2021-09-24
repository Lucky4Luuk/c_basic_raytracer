#ifndef _INCLUDE_DISPLAY_H_
#define _INCLUDE_DISPLAY_H_

#include <stdlib.h>
#include <stdio.h>

#if defined(_WIN32)
#include <Windows.h>
#elif defined(__linux__)
#include <sys/ioctl.h>
#endif // Windows/Linux

int SCREEN_WIDTH;
int SCREEN_HEIGHT;
uint8_t* SCREEN_BUF;

void get_terminal_size(int *width, int *height) {
#if defined(_WIN32)
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	int columns, rows;

	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
	columns = csbi.srWindow.Right - csbi.srWindow.Left + 1;
	rows = csbi.srWindow.Bottom - csbi.srWindow.Top + 1;
	*width = columns;
	*height = rows;
#elif defined(__linux__)
	struct winsize w;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
	*width = w.ws_col;
	*height = w.ws_row;
#endif // Windows/Linux
}

void createScreen() {
	int w = 0;
	int h = 0;

	get_terminal_size(&w, &h);

	// printf("cols: %d\n", cols);
	// printf("lines: %d\n", lines);

	SCREEN_WIDTH = w;
	SCREEN_HEIGHT = h;

	SCREEN_BUF = malloc(w * h * 3 * sizeof(*SCREEN_BUF));
	if (SCREEN_BUF == NULL) {
		fprintf(stderr, "Error allocating resources.\nProgram aborted.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT * 3; ++i) {
		SCREEN_BUF[i] = 0;
	}
}

void deleteScreen() {
	free(SCREEN_BUF);
}

void setPixel(int x, int y, uint8_t r, uint8_t g, uint8_t b) {
	int idx = (x + y * SCREEN_WIDTH) * 3;
	SCREEN_BUF[idx  ] = r;
	SCREEN_BUF[idx+1] = g;
	SCREEN_BUF[idx+2] = b;
}

void flushScreen() {
	for (int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT * 3; i+=3) {
		uint8_t r = SCREEN_BUF[i  ];
		uint8_t g = SCREEN_BUF[i+1];
		uint8_t b = SCREEN_BUF[i+2];
		printf("\x1b[38;2;%d;%d;%dm%c\x1b[0m", r,g,b, (char)219u);
	}
}

#endif
