/*! @file */
#ifndef HMPDF_OBJECT_H
#define HMPDF_OBJECT_H

/*! The data structure that all globally exposed functions receive
 *  a pointer to as their first argument.
 *  The user should not attempt to interact with this structure
 *  directly.
 */
typedef struct hmpdf_obj_s hmpdf_obj;

/*! Allocates a new hmpdf_obj.
 *  \return d       pointer to a new #hmpdf_obj,
 *                  NULL if memory allocation failed
 */
hmpdf_obj *hmpdf_new(void);

/*! Frees all memory associated with the hmpdf_obj.
 *  \param[in] d        #hmpdf_obj created with hmpdf_new()
 *  \return error code
 */
int hmpdf_delete(hmpdf_obj *d);

#endif
