#include "vector.h"
#include "common/assert.h"

#include <string.h>


void init(vector *v, size_t size_of_data) {
  v->size_of_data = size_of_data;
  v->size = 0;
  v->capacity = 1;
  v->data = malloc(v->size_of_data * v->capacity);
}

void *push(vector *v) {
  if (v->size == v->capacity) {
    v->capacity *= 2;
    v->data = realloc(v->data, v->size_of_data * v->capacity);
    RAY_VERIFY(v->data != NULL, "Failed to realloc");
  }

  v->size++;
  return back(v);
}

void *back(vector *v) {
  return (v->data + v->size_of_data * (v->size - 1));
}

