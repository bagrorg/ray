#pragma once

#include <stddef.h>

typedef struct {
  size_t size, capacity;
  size_t size_of_data;
  void *data;
} vector;

void init(vector *v, size_t size_of_data);
void *push(vector *v);
void *back(vector *v);


