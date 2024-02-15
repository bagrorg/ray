#ifndef ASSERT_H
#define ASSERT_H

#include <stdio.h>
#include <stdlib.h>

#define RAY_VERIFY(expr,...) \
	do { \
		if (!(expr)) { \
			fprintf(stderr, "Assertion failed at %s:%d: ", __FILE__, __LINE__);\
			fprintf(stderr, __VA_ARGS__); \
			fprintf(stderr, "\n"); \
			exit(EXIT_FAILURE); \
		}\
	} while (0)


#define RAY_WARNING(...) \
	do { \
		fprintf(stderr, "Warning at %s:%d: ", __FILE__, __LINE__);\
		fprintf(stderr, __VA_ARGS__); \
		fprintf(stderr, "\n"); \
		exit(EXIT_FAILURE); \
	} while (0)

#define RAY_FAIL(...) \
	do { \
		fprintf(stderr, "Execution failed at %s:%d: ", __FILE__, __LINE__);\
		fprintf(stderr, __VA_ARGS__); \
		fprintf(stderr, "\n"); \
		exit(EXIT_FAILURE); \
	} while (0)

#endif
