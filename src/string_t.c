#include <string.h>

#include "string_t.h"

/**
 * Initializes all fields of string_t struct. If str.buf is already
 * allocated, it is not freed.
 */
void string_set(string_t *str, const char *val) {
    size_t val_len = strlen(val);
    if (val_len >= string_t_initial_capacity) {
        str->buf = (char *) calloc(val_len*2, sizeof(char));
        str->capacity = val_len*2;
    } else {
        str->buf = (char *) calloc(string_t_initial_capacity, sizeof(char));
        str->capacity = string_t_initial_capacity;
    }
    str->size = val_len;
    memcpy(str->buf, val, val_len);
}

/**
 * Initializes a string_t structure. Writes to all fields.
 */
void string_init(string_t *str) {
    str->buf = (char *) calloc(string_t_initial_capacity, sizeof(char));
    str->capacity = string_t_initial_capacity;
    str->size = 0;
}

/**
 * Reserves the buffer of str to have provided capacity.
 *
 * If str already had a malloced buffer, it is not freed.
 */
void string_reserve(string_t *str, size_t capacity) {
    str->buf = (char *) calloc(capacity, sizeof(char));
    str->capacity = capacity;
    str->size = 0;
}


/**
 * Returns negative number on error.
 */
int string_appendc(string_t *str, char c) {
    if (str->size + 1 >= str->capacity) {
        size_t new_cap = 2 * str->capacity;
        if (new_cap < str->capacity) {
            // overflow
            return -2;
        }
        char *new_buf = (char *) realloc(str->buf, new_cap);
        if (new_buf == NULL) {
            return -1;
        }
        str->buf = new_buf;
        str->capacity = new_cap;
    }
    str->buf[str->size++] = c;
    str->buf[str->size] = 0x0;
    return 0;
}

/**
 * Returns negative number on error.
 */
int string_appends(string_t *str, const char *c) {
    int status;
    while ((*c) != 0x0) {
        status = string_appendc(str, *c);
        if (status < 0) {
            return status;
        }
        c++;
    }
    return 0;
}

// /**
//  * Returns negative number on error, 0 on success.
//  */
// int string_appends(string_t *str, string_t *ap) {

// }

/**
 * str.size should be checked before calling this. Calling
 * this on an empty string returns 0;
 */
char string_pop(string_t *str) {
    if (str->size > 0) {
        char c = str->buf[str->size-1];
        str->buf[str->size-1] = 0;
        str->size--;
        return c;
    }
    return 0;
}

/**
 * Sets the string size to 0.
 *
 * NOTE: does not erase all memory contents.
 */
void string_clear(string_t *str) {
    // memset(str->buf, 0, str->size);
    str->size = 0;
    str->buf[0] = 0x0;
}