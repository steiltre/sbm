#ifndef SBM_STACK_H
#define SBM_STACK_H

typedef struct sbm_ll_node {
  /* Value at node */
  int32_t val;

  /* Pointer to next node */
  struct sbm_ll_node * next;
} sbm_ll_node;

/*
 * @brief Create new node
 *
 * @return New node
 */
sbm_ll_node * sbm_ll_node_alloc();

/*
 * @brief Add new node to front of stack
 *
 * @param val Value node is to store
 * @param root Root node
 */
void sbm_ll_node_push(
    int32_t val,
    sbm_ll_node * root);

/*
 * @brief Remove top item from stack
 *
 * @param root Root node
 *
 * @return Value on top of stack
 */
int32_t sbm_ll_node_pop(
    sbm_ll_node * root);

#endif
