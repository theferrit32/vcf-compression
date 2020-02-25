#ifndef _STRING_T_H
#define _STRING_T_H

#include <stdlib.h>
#include <sys/types.h>

const size_t string_t_initial_capacity = 16;

typedef struct _string_t {
    char *buf;          // buffer
    size_t size;        // used space in the buffer
    size_t capacity;    // maximum size of the buffer;
} string_t;

void string_init(string_t *str);
void string_reserve(string_t *str, size_t capacity);
void string_set(string_t *str, const char *val);
int string_appendc(string_t *str, char c);
int string_appends(string_t *str, const char *c);
char string_pop(string_t *str);
void string_clear(string_t *str);

#endif