/*! @file hmpdf_data.h */
#ifndef HMPDF_DATA_H
#define HMPDF_DATA_H

/*! The data structure that all globally exposed functions receive
 *  a pointer to as their first argument.
 *  The user should not attempt to interact with this structure
 *  directly.
 */
typedef struct hmpdf_obj_s hmpdf_obj;

/*! Allocates a new hmpdf_obj.
 */
hmpdf_obj *hmpdf_new(void);

/*! Frees all memory associated with the hmpdf_obj.
 */
int hmpdf_delete(hmpdf_obj *d);

#endif
