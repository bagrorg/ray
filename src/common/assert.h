#ifndef ASSERT_H
#define ASSERT_H

#include <stdio.h>
#include <stdlib.h>

#define RAY_VERIFY(expr,...) \
	do { \
		if (!(expr)) { \
			fprintf(stderr, "Assertion failed: ");\
			fprintf(stderr, __VA_ARGS__); \
			fprintf(stderr, "\n"); \
			exit(EXIT_FAILURE); \
		}\
	} while (0)

#define RAY_FAIL(...) \
	do { \
		if (!(expr)) { \
			fprintf(stderr, "Execution failed: ");\
			fprintf(stderr, __VA_ARGS__); \
			fprintf(stderr, "\n"); \
			exit(EXIT_FAILURE); \
		}\
	} while (0)

#endif
